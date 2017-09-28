# BUGCOUNT README

----
## Usage
### Analysis by leaf insect count
    Rscript analysis.R WF_consolidated.tab \
                      leaf_resistance_analysis.pdf \
                      leaf;
                      
    Rscript blup_heritability.R WF_consolidated.tab\
            blup_analysis \
            leaf;

### Analysis by plant insect count
    Rscript analysis.R WF_consolidated.tab \
                      plant_resistance_analysis.pdf \
                      plant;
                      
    Rscript blup_heritability.R WF_consolidated.tab\
            blup_analysis \
            plant;
