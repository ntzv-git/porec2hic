# porec2hic

Découpe des reads Pore-C en monomères aux **vraies jonctions de ligation** (validées
par les reads HiFi), puis conversion **all-to-all** des monomères en paires pseudo-Hi-C
(R1/R2) pour les outils Hi-C (bwa mem -5SP, pairtools, YaHS, juicer…).

```bash
POREC_FILE=porec.fq.gz HIFI_FQ=hifi.fq.gz THREADS=96 ./porec2hic_hifi.sh
# -> porec2hic_out/porec_hic_R1.fastq.gz, porec_hic_R2.fastq.gz, porec_hic.stats.tsv
```

Dépendances : `minimap2`, `python3` (bibliothèque standard uniquement), `mawk`/`awk` ;
`pigz` est facultatif.

## Ce qui a changé par rapport au pipeline d'origine

### 1. Identification des jonctions avec les HiFi : alignement dans l'autre sens

Dans l'ancien pipeline, les reads HiFi (requêtes) étaient alignés sur les reads Pore-C
(cibles) en ne gardant que les alignements primaires (`-F 0x904`). Chaque read HiFi ne
couvre donc **qu'un seul** read Pore-C, alors qu'un locus est couvert par des dizaines
de reads Pore-C. La plupart des reads Pore-C se retrouvaient sans couverture HiFi : tous
leurs CATG passaient sous `MIN_COV`, et le pipeline faisait en réalité une digestion
in silico complète.

Mesure sur les données simulées de `tests/` (400 kb, HiFi 20x, 3 000 concatémères) :

| | reads Pore-C sans HiFi | précision des coupes | rappel |
|---|---|---|---|
| ancien pipeline | 72 % | **27 %** | 90 % |
| nouveau pipeline | 0,8 % | **96–98 %** | 89 % (94 % des jonctions dont le motif est intact) |

Le nouveau pipeline aligne les **Pore-C (requêtes) sur les HiFi (cibles)**. minimap2
coupe naturellement l'alignement à chaque jonction (primaire + supplémentaires), et
`-N` garde les autres reads HiFi du même locus (secondaires). Chaque site CATG est
ensuite classé à partir des intervalles alignés sur le read Pore-C :

| classe | règle | action |
|---|---|---|
| **P** protégé | ≥ `MIN_COV` alignements traversent `[x-FLANK, x+FLANK]` | CATG génomique, pas de coupe |
| **J** jonction | < `MIN_COV` traversants et ≥ `MIN_SIDE_COV` alignements d'un côté (fenêtre `SIDE_WINDOW`) | **coupe** |
| **j** voisin de jonction | autre CATG à moins de `2×FLANK` d'une jonction | pas de coupe, on garde le site le plus proche du point de cassure des alignements |
| **U** non résolu | pas d'information HiFi | pas de coupe (`CUT_UNRESOLVED=1` = ancien comportement) |

La différence importante avec `samtools depth` : la profondeur ne distingue pas un read
qui **traverse** le site d'un read qui **s'arrête** dessus, alors que les alignements
s'arrêtent justement aux jonctions. Le critère « traversant avec marge `FLANK` » fait
cette distinction.

`-N` est un plafond **global par read Pore-C**, pas par monomère : les secondaires du
monomère le plus long remplissent le quota. Avec `-N 10`, le rappel tombait à 47 %. La
valeur par défaut est maintenant 100 ; prévoir environ `MIN_COV × nb de monomères × 5`.

### 2. Coupe exactement au site de restriction

- `CUT_OFFSET` est la position de coupure **dans** le motif : NlaIII `CATG^` = 4,
  DpnII `^GATC` = 0. L'ancien `bedtools subtract` supprimait une base (le G) à chaque
  coupure.
- `DUP_MOTIF=1` (par défaut) : le CATG, reconstitué à la ligation, est conservé sur les
  **deux** monomères (`…CATG | CATG…`). Chaque monomère correspond alors exactement à
  sa séquence génomique, quelle que soit l'orientation du fragment ligué.
- Motifs IUPAC acceptés. Les motifs non palindromiques sont recherchés sur les deux brins.
- Les sites à moins de `FLANK` pb d'une extrémité de read ne sont pas coupés, et les
  monomères de moins de `MIN_MONO_LEN` (50 pb) sont écartés car non mappables.

### 3. Regroupement all-to-all en une seule passe

Les étapes 2 à 7 (seqkit locate, samtools depth, bedtools, seqkit subseq, awk) sont
remplacées par un seul passage en flux (`porec_hifi_split.py`) :

- pas de BAM, pas de tri ni d'index, pas de fichiers BED/depth intermédiaires de
  plusieurs milliards de lignes, pas de passage `seqkit stats`/`fx2tab` ;
- le PAF est compacté à la volée en une ligne par read (`read  qlen  starts  ends`) ;
- `THREADS-1` workers ; chaque worker écrit ses propres fragments gzip, concaténés à la
  fin (une concaténation de membres gzip reste un gzip valide). R1 et R2 restent
  synchronisés et aucune séquence ne repasse par l'IPC ;
- pour un read à n monomères, toutes les paires C(n,2) sont écrites, nommées
  `@read:i-j/1` et `@read:i-j/2` ;
- `MAX_MONOMERS` permet d'écarter les concatémères extrêmes (n² paires).

## Paramètres principaux

| variable | défaut | rôle |
|---|---|---|
| `MOTIF` / `CUT_OFFSET` | `CATG` / `4` | enzyme |
| `MIN_COV` | 3 | alignements traversants pour protéger un site |
| `MIN_SIDE_COV` | `MIN_COV` | alignements d'un côté pour valider une jonction |
| `FLANK` / `SIDE_WINDOW` | 25 / 100 | marges (pb) |
| `MM2_PRESET` | `map-ont` | `lr:hq` pour ONT R10 Q20+ (minimap2 ≥ 2.27) |
| `MM2_N` | 100 | secondaires minimap2 |
| `MM2_BATCH` | 50G | taille des lots d'index HiFi (`-I`, avec `--split-prefix`) : environ 2–2,5 octets/base de RAM |
| `MM2_EXTRA` | – | options minimap2 supplémentaires (ex. `-c` pour des bornes exactes, plus lent) |
| `WRITE_SITES=1` | 0 | table par site (`span`, `covL`, `covR`, classe) pour calibrer `MIN_COV` |
| `WRITE_MONOMERS=1` | 0 | FASTQ des monomères |

Chaque étape écrit un fichier `stepN.done` : une relance reprend là où le pipeline s'est
arrêté.

## Remarques

- **Temps et mémoire de l'étape 1** : l'index est construit sur les HiFi. Pour
  `MIN_COV=3`, une couverture HiFi de 10–15x suffit. Sous-échantillonner les HiFi
  (`seqkit sample`) réduit fortement la mémoire et le temps.
- Les jonctions entre fragments génomiquement proches (< ~5 kb, même orientation)
  peuvent être chaînées par minimap2 à travers un seul read HiFi, et restent alors
  protégées. Ce sont des contacts cis très courts, peu informatifs.
- Une jonction dont le motif a été altéré par une erreur de séquençage n'est pas coupée,
  puisque la coupe se fait au motif.

## Test

```bash
pip install mappy
tests/run_test.sh /tmp/porec2hic_test   # simulation + pipeline + précision/rappel
```
