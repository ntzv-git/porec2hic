# porec2hic

Découpe des reads Pore-C en monomères aux **vraies jonctions de ligation**, validées par
les reads HiFi, puis conversion **all-to-all** des monomères en paires pseudo-Hi-C
(R1/R2).

```bash
POREC_FILE=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
# -> porec2hic_out/porec_hic_R1.fastq.gz, porec_hic_R2.fastq.gz, sites.final.tsv.gz
```

Outils : `minimap2`, `seqkit`, `bedtools` (≥ 2.26), `gawk`, `pigz` (facultatif).
Le pipeline ne contient aucun script Python.

## Étapes

| étape | outil | sortie |
|---|---|---|
| 1. longueurs des reads | `seqkit fx2tab -n -l` | `porec.genome` (ordre du FASTQ) |
| 2. alignement Pore-C → HiFi | `minimap2 -c -N` + `gawk` | `hifi_blocks.bed.gz` : blocs continus de chaque read HiFi, en coordonnées Pore-C |
| 3. sites de restriction | `seqkit locate` + `awk` | position de coupure exacte `x = début motif + CUT_OFFSET` |
| 4. continuité HiFi au site | `bedtools map` (×3) | `span`, `covL`, `covR` |
| 5. un seul site par jonction | `bedtools map -o collapse` + `gawk` | `sites.final.tsv.gz` |
| 6. découpe + paires all-to-all | `seqkit fx2tab` + `awk` | R1/R2 |

Tous les fichiers restent dans l'ordre du FASTQ Pore-C, et `bedtools map -g porec.genome`
travaille dans cet ordre : il n'y a aucun tri global. Chaque étape écrit un fichier
`stepN.done`, donc une relance reprend là où le pipeline s'est arrêté.

## Définition d'un read HiFi « continu » à travers un site

Les reads Pore-C (requêtes) sont alignés sur les reads HiFi (cibles), avec `-N` pour
garder les autres reads HiFi du même locus. Chaque alignement est ensuite converti en
**blocs continus**, en coordonnées sur le read Pore-C :

1. l'alignement est coupé à chaque indel ≥ `MAX_GAP` (lu dans le CIGAR, option `-c`) ;
2. les blocs d'un **même read HiFi**, sur le même brin, qui sont colinéaires et séparés
   de moins de `MAX_GAP` sur le Pore-C **et** sur le HiFi sont fusionnés.

Conséquences pour un read HiFi aligné sur les bases qui flanquent un site :

| situation | lecture | effet |
|---|---|---|
| un seul bloc couvre `[x-FLANK, x+FLANK]` | le HiFi lit la séquence à travers le site : CATG génomique | compte pour la protection |
| deux blocs du même read HiFi, colinéaires, trou < `MAX_GAP` (erreur de séquençage Pore-C sur le site) | fusionnés : continu | compte pour la protection |
| deux blocs du même read HiFi avec un saut ≥ `MAX_GAP`, ou sur des brins différents (ligation cis à courte distance, inversion) | non continu | compte seulement comme couverture latérale |
| blocs de reads HiFi différents, l'un finit avant `x` et l'autre commence après | signature d'une jonction | couverture latérale (`covL`/`covR`) |
| un bloc dépasse `x` de moins de `FLANK` pb | débordement d'alignement, pas une preuve de continuité | couverture latérale |

## Classes des sites

| classe | condition | action |
|---|---|---|
| **P** protégé | `span ≥ MIN_COV` : au moins `MIN_COV` reads HiFi **distincts** ont un bloc continu couvrant entièrement `[x-FLANK, x+FLANK]` | pas de coupe |
| **J** jonction | `span < MIN_COV` et au moins `MIN_SIDE_COV` reads HiFi distincts ont un bloc dans `[x-SIDE_WINDOW, x)` ou dans `[x, x+SIDE_WINDOW)` | **coupe** |
| **j** motif voisin d'une jonction | site J à moins de `2×FLANK` d'un autre site J : on garde celui le plus proche du point de cassure des blocs HiFi | pas de coupe |
| **U** non résolu | pas assez de reads HiFi autour du site | pas de coupe (`CUT_UNRESOLVED=1` pour couper) |

Aux extrémités d'un read Pore-C, la fenêtre est tronquée à la séquence disponible.

## Choix de `FLANK`

À une jonction NlaIII, le motif CATG appartient aux **deux** fragments ligués. Un
alignement HiFi dépasse donc la jonction d'au moins 4 pb, plus quelques bases qui
correspondent par hasard. `FLANK` doit être supérieur à ce débordement.

Débordement des blocs HiFi au-delà d'une vraie jonction (simulation, minimap2 2.28 `-c`) :

| débordement (pb) | 1 | 2 | 3 | **4** | 5 | 6 | 7 | 8 | 9 | 10 | 11–15 | > 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| % des blocs | 14 | 5 | 2 | **56** | 14 | 4 | 1,2 | 1,2 | 0,8 | 0,5 | 1 | < 0,5 |

Résultat selon `FLANK` (précision et rappel des coupes) :

| FLANK | 0 | 5 | 10 | **15** | 20 | 25 | 50 |
|---|---|---|---|---|---|---|---|
| précision | 71 % | 97,4 % | 97,9 % | **97,9 %** | 97,8 % | 97,7 % | 97,2 % |
| rappel | 5 % | 69 % | 91,8 % | **92,7 %** | 92,1 % | 91,4 % | 87,6 % |

- Une fenêtre trop petite fait protéger les vraies jonctions (débordement).
- Une fenêtre trop grande fait tomber davantage de motifs voisins d'une jonction sous le
  seuil ; la plupart sont rattrapés à l'étape 5 (classe `j`).
- Valeur par défaut : 15 pb. Le rappel manquant correspond surtout aux jonctions dont le
  motif a été altéré par une erreur de séquençage (≈ 4 %), qui ne peuvent pas être
  coupées au motif.

## Découpe et paires

- `CUT_OFFSET` est la position de coupure dans le motif : NlaIII `CATG^` = 4,
  DpnII `^GATC` = 0.
- `DUP_MOTIF=1` : le motif reconstitué à la ligation est gardé sur les deux monomères
  (`…CATG | CATG…`).
- **Aucun filtrage** : tous les monomères sont gardés, quelle que soit leur taille, et
  un read de n monomères donne n(n-1)/2 paires `@read:i-j/1` et `@read:i-j/2`. Les
  filtres (taille, MAPQ, distance…) se font en aval sur les paires.
- Avec `WRITE_MONOMERS=1`, le FASTQ des monomères est aussi écrit.

## Paramètres

| variable | défaut | rôle |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / 4 | enzyme (motif palindromique, IUPAC accepté) |
| `MIN_COV` | 3 | reads HiFi continus nécessaires pour protéger un site |
| `FLANK` | 15 | marge de continuité (pb) |
| `MAX_GAP` | 50 | indel/trou (pb) qui interrompt la continuité |
| `MIN_SIDE_COV` / `SIDE_WINDOW` | `MIN_COV` / 100 | preuve HiFi à côté d'une jonction |
| `CUT_UNRESOLVED` | 0 | couper aussi les sites U |
| `MM2_PRESET` | `map-ont` | `lr:hq` pour ONT R10 Q20+ |
| `MM2_N` | 100 | secondaires minimap2 ; plafond **par read Pore-C**, prévoir ≈ `MIN_COV` × nb de monomères × 5 |
| `MM2_BATCH` | 50G | lots d'index HiFi (`-I` avec `--split-prefix`) |

## Coût

Le volume d'alignements est d'environ (nombre de monomères Pore-C) × (profondeur HiFi),
avec CIGAR. Pour `MIN_COV=3`, 10–15x de HiFi suffisent : sous-échantillonner
(`seqkit sample`) réduit d'autant le temps de minimap2 et la taille de
`hifi_blocks.bed.gz`.

## Test

```bash
tests/run_test.sh /tmp/porec2hic_test   # simulation + pipeline + précision/rappel
```

La simulation et l'évaluation (`tests/*.py`) sont en Python, mais elles ne font pas
partie du pipeline.
