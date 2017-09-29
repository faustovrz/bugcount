# BUGCOUNT README

----
## Usage
### Analysis by plant insect count
    Rscript analysis.R  plant resistance_analysis wf_consolidated.tab \
            > plant_resistance_analysis.log 2>&1
    
    Rscript blup_heritability.R plant blup_analysis wf_consolidated.tab \
            > plant_blup_analysis.log 2>&1

### Analysis by leaf insect count
    Rscript analysis.R leaf resistance_analysis wf_consolidated.tab \
            > leaf_resistance_analysis.log 2>&1
            
    Rscript blup_heritability.R leaf blup_analysis wf_consolidated.tab \
            > leaf_blup_analysis.log 2>&1

