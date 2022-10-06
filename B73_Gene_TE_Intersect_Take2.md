
### Breakdown of Strategy: 
1. a curated list of B73 structural TEs was used (no helitrons) 
2. the canonical transcripts from B73 were used to generate gene and intron files 
3. I determined any TE that overlaps coding sequence and deleted it for downstream analysis 
4. Remaining TEs were intersected with Intron locations to generate a list of TEs that overlap Introns 
    Note that there are introns that overlap both coding and intron - these are classed as coding here 
5. Remaining TEs were seperately intersected with 2kb upstream (promoter) and 2kb down files 
    I chose 2kb this time based on a recent Brassica paper https://onlinelibrary.wiley.com/doi/10.1111/pbi.13807
    
    This Brassica paper has the following figure - I am recreating this for B73\

<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184729988-cfd0d635-2884-4654-b404-335fb267a11d.png?raw=true" alt="image"/>
</p>

### Modifications to Strategy:
1. I used Manisha's curated structural and homology TE list instead of Structural only (though I did a grouped and seperate analysis)
2. Canonical transcripts were still used to generate the Intron files
3. I did Intron specific overlaps, however for the specific below analysis I'm just saying "does a TE overlap a canonical gene"
4. I created files with 500bp blocks up and downstream of each gene - 4 blocks each, so looking at 2kb up and down \
Note that I did not correct for any genes that are within 4kb of eachother, so there may be some genomic regions that get counted more than once.

To generate the below figure I calculated the number of each class of TE 

![image](https://user-images.githubusercontent.com/43852873/194428494-63168df5-4f10-4b75-9ffb-f08bdf72cde3.png)


### Here is the B73 data

### ALL GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731583-ab2c567e-7b18-4ba6-a4eb-60b6953b7f9b.png?raw=true" alt="image"/>
</p>

### NLR GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731759-578c3e9c-f3d4-448f-a730-2b7cb427bfcd.png?raw=true" alt="image"/>
</p>


## Reorganizing the data
This is interesting but it's challenging because the datasets are such different sizes - Manisha suggested pulling 2 more subsets of genes \
and seeing what their patterns look like

<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/185230446-e4eaadb2-e142-4369-9f0a-9cd997697f2c.png?raw=true" alt="image"/>
</p>


