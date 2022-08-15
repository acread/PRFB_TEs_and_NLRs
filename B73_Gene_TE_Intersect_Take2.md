
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

Lots of code....

### Here is the B73 data

### ALL GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731583-ab2c567e-7b18-4ba6-a4eb-60b6953b7f9b.png?raw=true" alt="image"/>
</p>

### NLR GENES
<p align="center">
  <img src="https://user-images.githubusercontent.com/43852873/184731759-578c3e9c-f3d4-448f-a730-2b7cb427bfcd.png?raw=true" alt="image"/>
</p>

