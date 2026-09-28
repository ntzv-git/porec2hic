# porec2hic

Découpe des reads Pore-C en monomères aux **vraies jonctions de ligation**, validées par
les reads HiFi, puis conversion **all-to-all** des monomères en paires pseudo-Hi-C
(R1/R2).

```bash
POREC_FQ=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
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
| 4. protection HiFi | `bedtools map -f 1.0 -o count_distinct` | `sites.final.tsv.gz` (`span`, classe P/C) |
| 5. découpe, filtre des pseudo-monomères, paires all-to-all | `seqkit fx2tab` + `awk` | R1/R2 |

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
| deux blocs du même read HiFi avec un saut ≥ `MAX_GAP`, ou sur des brins différents (ligation cis à courte distance, inversion) | non continu | ne protège pas |
| blocs de reads HiFi différents, l'un finit avant `x` et l'autre commence après | signature d'une jonction | ne protège pas |
| un bloc dépasse `x` de moins de `FLANK` pb | débordement d'alignement, pas une preuve de continuité | ne protège pas |

## Règle de coupure

| classe | condition | action |
|---|---|---|
| **P** protégé | `span ≥ MIN_COV` : au moins `MIN_COV` (3) reads HiFi **distincts** ont un bloc continu couvrant entièrement `[x-FLANK, x+FLANK]` | pas de coupe |
| **C** coupé | `span < MIN_COV`, y compris quand aucun read HiFi n'est présent | **coupe** |

Aux extrémités d'un read Pore-C, la fenêtre est tronquée à la séquence disponible.

Cette règle privilégie les coupures : un site mal couvert par les HiFi est coupé. Sur la
simulation, 99,4 % des jonctions dont le motif est intact sont coupées (rappel 95,7 %
en comptant les motifs altérés par une erreur de séquençage). Il y a 787 coupes sur des
CATG génomiques, pour 7 295 vraies jonctions coupées :

- 552 sont à moins de 20 pb d'une vraie jonction (voir ci-dessous) ;
- 235 sont ailleurs, sur des sites couverts par 0 à 2 reads HiFi seulement.

## Pseudo-monomères

Un CATG situé à moins de ~`FLANK` pb d'une vraie jonction n'est traversé par aucun read
HiFi avec la marge requise. Il est donc coupé lui aussi, et entre les deux coupes naît
un **pseudo-monomère** de quelques pb. Taille des monomères sur la simulation :

| longueur (pb) | 0–4 | 5–9 | 10–14 | 15–19 | 20–24 | 25–29 | 30–49 |
|---|---|---|---|---|---|---|---|
| pseudo-monomères (artefacts) | 4 | 132 | 304 | 131 | 18 | 16 | 69 |
| vrais monomères | 0 | 24 | 63 | 76 | 43 | 35 | 185 |

Leur taille maximale vaut environ `FLANK` + longueur du motif (19 pb pour CATG). Les
monomères de moins de `MIN_MONO_LEN` (par défaut `FLANK` + longueur du motif + 1 = 20 pb)
sont donc écartés **avant** l'appariement. Les vrais monomères de cette taille ne sont de
toute façon pas mappables de manière unique. Pour ne rien filtrer, utiliser
`MIN_MONO_LEN=1`.

## Choix de `FLANK`

À une jonction NlaIII, le motif CATG appartient aux **deux** fragments ligués. Un
alignement HiFi dépasse donc la jonction d'au moins 4 pb, plus quelques bases qui
correspondent par hasard. `FLANK` doit être supérieur à ce débordement.

Débordement des blocs HiFi au-delà d'une vraie jonction (simulation, minimap2 2.28 `-c`) :

| débordement (pb) | 1 | 2 | 3 | **4** | 5 | 6 | 7 | 8 | 9 | 10 | 11–15 | > 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| % des blocs | 14 | 5 | 2 | **56** | 14 | 4 | 1,2 | 1,2 | 0,8 | 0,5 | 1 | < 0,5 |

Résultat selon `FLANK` (précision et rappel des coupes, mesurés avec la version
précédente de la règle de coupure) :

| FLANK | 0 | 5 | 10 | **15** | 20 | 25 | 50 |
|---|---|---|---|---|---|---|---|
| précision | 71 % | 97,4 % | 97,9 % | **97,9 %** | 97,8 % | 97,7 % | 97,2 % |
| rappel | 5 % | 69 % | 91,8 % | **92,7 %** | 92,1 % | 91,4 % | 87,6 % |

- Une fenêtre trop petite fait protéger les vraies jonctions (débordement).
- Une fenêtre trop grande fait tomber davantage de motifs voisins d'une jonction sous le
  seuil, ce qui crée des pseudo-monomères plus longs.
- Valeur par défaut : 15 pb. Le rappel manquant correspond surtout aux jonctions dont le
  motif a été altéré par une erreur de séquençage (≈ 4 %), qui ne peuvent pas être
  coupées au motif.

## Découpe et paires

- `CUT_OFFSET` est la position de coupure dans le motif : NlaIII `CATG^` = 4,
  DpnII `^GATC` = 0.
- `DUP_MOTIF=1` : le motif reconstitué à la ligation est gardé sur les deux monomères
  (`…CATG | CATG…`).
- Seul filtrage : les pseudo-monomères de moins de `MIN_MONO_LEN`. Tous les autres
  monomères sont gardés, et un read de n monomères retenus donne n(n-1)/2 paires
  `@read:i-j/1` et `@read:i-j/2`. Les autres filtres (MAPQ, distance…) se font en aval
  sur les paires.
- Avec `WRITE_MONOMERS=1`, le FASTQ des monomères retenus est aussi écrit.

## Paramètres

| variable | défaut | rôle |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / 4 | enzyme (motif palindromique, IUPAC accepté) |
| `MIN_COV` | 3 | reads HiFi continus nécessaires pour protéger un site ; sinon coupure |
| `FLANK` | 15 | marge de continuité (pb) |
| `MAX_GAP` | 50 | indel/trou (pb) qui interrompt la continuité |
| `MIN_MONO_LEN` | `FLANK` + motif + 1 (20) | monomères écartés avant l'appariement |
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
