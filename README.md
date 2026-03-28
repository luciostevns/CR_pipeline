**Computationally predicting T-cell cross-reactivity**

**induced autoimmunity with pathogenic proteome**



The purpose of the following project is to locate possible sequence

alignment between pathogenic proteins and known autoimmune epitopes.





DATA COLLECTION:



* Collect data from IEDB.org with epitopes and metadata of interest, this should at the vary least contain protein sequences.

   Current script is specific for the following metadata: (Assay ID - IEDB IRI, Epitope - Name, Epitope - Molecule Parent,

   Epitope - Molecule Parent IRI, 1st in vivo Process - Disease, 1st in vivo Process - Disease Stage, MHC Restriction - Name,

   Epitope - Starting Position, Epitope - Ending Position Epitope - Modified residues). If not change script to fit metadata



* Collect data from uniport.org with proteomes of pathogens of interest and download .tsv file.

   Including Assembly id is important.

 

* Collect data from ncbi.nlm.nih.gov/pathogens/isolates (This is for verifying chosen pathogens are pathogenic).

   Including Assembly id is important.







CODE STUFF (REPLACE WITH BETTER HEADER)



1\. Clean data and filter (./R\_stuff/IEDB\_wrangling.R). This will filter protein of being in

   specific length range (12-25), it will remove if one epitope is nested in a bigger

   epitope, making it redundant. Then it will generate all possible 9mers of each epitope

   it will then write a .csv file, with chosen columns as columns and each row being a 9mer





 

