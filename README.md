<h2>Identifiability of word models</h2>



To run the INDRA assembly code, use the following command,in case the policy type is passed as a dictionary i.e. different and specific policies for different events
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy '{\"Phosphorylation\":\"two_step\",\"Dephosphorylation\":\"two_step\"}'
```
Use the followint command, in case the policy type is global and has to be passed as a single string. For example "one_step", "two_step", "michaelis_menten"
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy "one_step"
```

We have considered five cases on which we have applied assembling of INDRA model using REACH and compared the results with the TRIPS.
The input by the user has to be of the format --idtypes "IDTYPE1" "IDTYPE2"


Two output files will be saved in the output folder. output.txt specifies the desired dynamical model in the DAISY format while attributes.txt mentions the attributes of the assembled model like statements, rules, monomers and the number of ODEs in the reaction system. In case you want to change the name of the output files, you can specify them as follows
```
ython main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy "one_step" --output_txt "new_output_name.txt" --attr_txt "new_attribute_name.txt"
```
