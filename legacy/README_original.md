# bacdive_harvester

An R toolkit for extracting and analyzing bacterial trait data from the BacDive database (Bacterial Diversity Metadatabase). This tool automates the collection of morphological, physiological, and metabolic characteristics from bacterial records.

For an input list of species names this toolkit can extract:
- Morphological traits: Gram stain, cell shape, motility
- Metabolic profiles: Positive metabolite utilization, enzyme activities
- Physiological characteristics: Oxygen tolerance (aerobic/anaerobic classification)
- Strain information: Type strain identification and strain designations

In the example provided it pulls species names from the Human Microbiome Project via the curatedMetagenomicData BioConductor package. 

BacDive API Access
⚠️ Important: You must sign up for free API access to use this tool.

- Visit https://api.bacdive.dsmz.de
- Create a free account
- Obtain your API credentials (username and password)
- Add your credentials to the script:
