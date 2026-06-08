
# Superimposition Tutorial

## Introduction
This tutorial outlines the steps of performing systematic superimosition of proteins. The process of superimposition is useful in comparing structures, analyzing active sites and their D structures.

The `superimposition.py` script accounts for 3 scenarios.
+ Same proteins but different variant with same number of atoms
+ Proteins of different lengths either cause they are diffent proteins or a protein mutant with insertions/deletions or missing residues
+ Same length but different number of atoms

In the event that the proteins are different lengths the script defaults to performing alignment on the *alpha* carbons. CA alignment however, is also an option the user can specify when writing the command in the terminal. For the last scenario the user can skip alignment and only perform superimposition. Finally, using an CLI flag the user can visualize the pairwise alignment with match lines.

For the tutorial we compare the enzymes and variants of FLS (formolase), and BAL (benzaldehyde lyase). Specifically PfBAL from _Pseudomonas fluorescens_ to demonstrate the three scenarios. FLS is a 7-site mutant of PfBAL so we can expect the RMSD scores between the proteins to be of interest.


## Environment setup
The main package required for the tutorial is Biopython. With the command below you can create an environment dedicated to superimposition, it includes Biopython and a couple other useful packages.

```
conda env create --file superimpose_env.yml
```

Activate your environment:
```
conda activate superimpose
```
## Script
Here are all the options the script supports and how to see them:
```
python superimposition.py --help
```
Output:
```
usage: superimposition.py [-h] [--chain1 CHAIN1] [--chain2 CHAIN2] [--ca] [-o OUTPUT] [--show-alignment-v] [--no-align] pdb1 pdb2

Superimpose two protein structures using Biopython

positional arguments:
  pdb1                  Reference PDB
  pdb2                  Mobile PDB

options:
  -h, --help            show this help message and exit
  --chain1 CHAIN1       Chain in reference structure
  --chain2 CHAIN2       Chain in mobile structure
  --ca                  Force alpha carbons only
  -o OUTPUT, --output OUTPUT
                        Write transformed PDB
  --show-alignment-v    Print the sequence alignment (target/query with match line)
  --no-align            Skip sequence alignment and superimpose residues by index
```

The only required arguments are the structre PDB files, all else is optional.
As shown in the help message the first structure input, is the one used as **reference**. 


## Structures
### 1. Same protein, same length
For this example we will use FLS `4QQ8.pdb`, since it is a homotetramer we will use our script to align and superimpose two chains of the enzyme.

First to get the sequence:

```
wget https://files.rcsb.org/download/4QQ8.pdb
```

+ To align the **chain structures**:
```
python superimposition.py 4QQ8.pdb 4QQ8.pdb --chain1 A --chain2 B 
```

**Expected output:**
```
Aligned atoms: 4152
RMSD: 0.403 Å
```

+ Or to align only based on the **_alpha_ carbons** of each residue when adding the --ca flag:
```
python superimposition.py 4QQ8.pdb 4QQ8.pdb --chain1 A --chain2 B --ca -o 4qq8-chainA-chainB-output.pdb
```
**Expected output:**
```
Aligned atoms: 563
RMSD: 0.194 Å
Saved superimposed structure: 4qq8-chainA-chainB-output.pdb
```
Here we have chosen to also add an output PDB file. The **output PDB** contains the transformed mobile PDB, that is the mobile sequence rotated and translated into the reference coordinate frame.

---

### 2. Proteins of different lenghts
For this example we will once again use FLS `4QQ8.pdb` and a PfBAL monomer `pdb2ag0.ent` (for the specific PfBAL variant we will see in the output that it is 554 residues long, the sequence is actually 563 residues long however in EMBL-EBI the last few residues are listed as unobserved sequences). \
For PfBAL we will use a variant of [Q9F4L3](https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/Q9F4L3). The structure can be downloaded is a PBD file with an .ent extension. [The .ent extension does not affect the file or its processing here]

+ For Q9F4L3-2ag0
```
wget https://www.ebi.ac.uk/pdbe/entry-files/download/pdb2ag0.ent
```

To superimpose the two structures:

```
python superimposition.py pdb2ag0.ent 4QQ8.pdb
```

**Expected output:**
```
Sequences have different lengths:
  Structure 1: 554 residues
  Structure 2: 563 residues
Switching automatically to CA-only superposition.

Aligned atoms: 554
RMSD: 0.566 Å
```
+ Here it might also be helpful to visialize the alignment.\
 For that we can add the `--show-alignment-v` flag:
```
python superimposition.py pdb2ag0.ent 4QQ8.pdb --show-alignment-v
```

**Expected output:**
```
  Structure 1: 554 residues
  Structure 2: 563 residues
Switching automatically to CA-only superposition.


Sequence alignment:

target            0 AMITGGELVVRTLIKAGVEHLFGLHGAHIDTIFQACLDHDVPIIDTRHEAAAGHAAEGYA
                  0 ||||||||||||||||||||||||||.|||||||||||||||||||||||||||||||||
query             0 AMITGGELVVRTLIKAGVEHLFGLHGIHIDTIFQACLDHDVPIIDTRHEAAAGHAAEGYA

target           60 RAGAKLGVALVTAGGGFTNAVTPIANAWLDRTPVLFLTGSGALRDDETNTLQAGIDQVAM
                 60 |||||||||||||||||||||||||||..|||||||||||||||||||||||||||||||
query            60 RAGAKLGVALVTAGGGFTNAVTPIANARTDRTPVLFLTGSGALRDDETNTLQAGIDQVAM

target          120 AAPITKWAHRVMATEHIPRLVMQAIRAALSAPRGPVLLDLPWDILMNQIDEDSVIIPDLV
                120 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query           120 AAPITKWAHRVMATEHIPRLVMQAIRAALSAPRGPVLLDLPWDILMNQIDEDSVIIPDLV

target          180 LSAHGARPDPADLDQALALLRKAERPVIVLGSEASRTARKTALSAFVAATGVPVFADYEG
                180 ||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||
query           180 LSAHGAHPDPADLDQALALLRKAERPVIVLGSEASRTARKTALSAFVAATGVPVFADYEG

target          240 LSMLSGLPDAMRGGLVQNLYSFAKADAAPDLVLMLGARFGLNTGHGSGQLIPHSAQVIQV
                240 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query           240 LSMLSGLPDAMRGGLVQNLYSFAKADAAPDLVLMLGARFGLNTGHGSGQLIPHSAQVIQV

target          300 DPDACELGRLQGIALGIVADVGGTIEALAQATAQDAAWPDRGDWCAKVTDLAQERYASIA
                300 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query           300 DPDACELGRLQGIALGIVADVGGTIEALAQATAQDAAWPDRGDWCAKVTDLAQERYASIA

target          360 AKSSSEHALHPFHASQVIAKHVDAGVTVVADGALTYLWLSEVMSRVKPGGFLCHGYLGSM
                360 ||||||||||||||||||||||||||||||||.||||||||||||||||||||||||.||
query           360 AKSSSEHALHPFHASQVIAKHVDAGVTVVADGGLTYLWLSEVMSRVKPGGFLCHGYLNSM

target          420 GVGFGTALGAQVADLEAGRRTILVTGDGSVGYSIGEFDTLVRKQLPLIVIIMNNQSWGAT
                420 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|
query           420 GVGFGTALGAQVADLEAGRRTILVTGDGSVGYSIGEFDTLVRKQLPLIVIIMNNQSWGWT

target          480 LHFQQLAVGPNRVTGTRLENGSYHGVAAAFGADGYHVDSVESFSAALAQALAHNRPACIN
                480 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
query           480 LHFQQLAVGPNRVTGTRLENGSYHGVAAAFGADGYHVDSVESFSAALAQALAHNRPACIN

target          540 VAVALDPIPPEELI--------- 554
                540 ||||||||||||||--------- 563
query           540 VAVALDPIPPEELILIGMDPFAG 563

Aligned atoms: 554
RMSD: 0.566 Å
```
>This output however might be too long and inconvenient to have in the terminal. Instead you can redirect the output to a separate .txt file:
```
python superimposition.py pdb2ag0.ent 4QQ8.pdb --show-alignment-v > testoutput.txt
```

---
### 3. Same length - different number of atoms - no alignment
Here we introduce another variant of PfBAL, [P51853](https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/P51853) to compare to [Q9F4L3](https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/Q9F4L3). 

For P51853-2uz1
```
wget https://www.ebi.ac.uk/pdbe/entry-files/download/pdb2uz1.ent
```


We can expect that even though these sequences are the same length they will have a different number of atoms.
+ We can choose to skip alignment and only superimpose:

```
python superimposition.py pdb2uz1.ent pdb2ag0.ent --no-align
```
**Expected output:**
```
Skipping sequence alignment. Superimposing residues by index.

554 residues superimposed
RMSD: 0.806 Å
```


+ Or we can still choose to forego alignment but still match based on the _alpha_ carbons:

```
python superimposition.py pdb2uz1.ent pdb2ag0.ent --ca --no-align
```
**Expected output:**
```
Skipping sequence alignment. Superimposing residues by index.

Aligned atoms: 554
RMSD: 0.660 Å
```

---

### Result Interpretation
| RMSD  | Meaning                     |
| ----- | --------------------------- |
| < 1 Å | nearly identical structures |
| 1–2 Å | same fold                   |
| 2–4 Å | homologous                  |


+ You can compare the resulting PDB up against the PDB of the reference structure and look at columns 7-9 for the x, y and z coordinates of the atom. The less these values deviate from pdb to pdb the better aligned the two structures are.

---

To visualize the superimposition using *PyMOL*:
```
load 4QQ8.pdb, ref
load 4qq8-chainA-chainB-output.pdb, mob
```

## Future directions
+ Base the RMSD on region sequence alignment

+ In the case of multimers the user could choose to superimpose selected chains

+ Include outlier rejection 

+ Setting of a default aligned atom Å value, but include the possibility of specifying a specific Å cutoff value

+ The script could eventually expand to use as input a refernce sequence and a list of proteins to compare it to, where the output will be a csv file compliling the RMSD scores

## Resources

[Biopython](https://biopython.org/)
