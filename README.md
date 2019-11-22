# CodingFrame
A tool to identify which of the six frames of a DNA backbone is the protein coding frame


The current code reads in a fasta formatted genome file, efficiently parses it to calculate the
entropy of three forward frames in amino-acid space.  It them uses the matplotlib library to plot
the estimated entropy at each base position along the genome, of the the three forward frames.

## To Run
```sh
python3 codingframe.py FASTA_FILE
```


## Testing
For testing purposes we choose the genome 1007869.8, due to it's smaller size (~140 kilo bases),
and that the orientation of the fasta file puts it's genes predominantly in the forward direction. 

```sh
python3 codingframe.py genomes/fna/1007869.8.fna.gz
```

After running the above command you should get an output image similiar to the one below
![](https://github.com/deprekate/CodingFrame/blob/master/figure.png)
