## Setup

```bash
genhub-build.py --cfgdir=config/ \
                --workdir=species/ \
                --numprocs=10 \
                --batch=modorg \
                download format prepare stats
```

## Tables

```bash
for species in Scer Cele Crei Mtru Dmel Agam Xtro Drer Mmus Hsap
do
    python pilocus_summary.py --gff3=species/${species}/${species}.iloci.gff3 \
                              species/${species}/${species}.iloci.tsv \
                              species/${species}/${species}.pre-mrnas.tsv \
        | tail -n 1
done

for species in Scer Cele Crei Mtru Dmel Agam Xtro Drer Mmus Hsap
do
    python nilocus_summary.py --gff3=species/${species}/${species}.iloci.gff3 \
                              species/${species}/${species}.iloci.tsv \
        | tail -n 1
done

for species in Scer Cele Crei Mtru Dmel Agam Xtro Drer Mmus Hsap
do
    python iilocus_summary.py --gff3=species/${species}/${species}.iloci.gff3 \
                              species/${species}/${species}.iloci.tsv \
        | tail -n 1
done
```
