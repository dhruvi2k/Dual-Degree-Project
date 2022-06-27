<h2 align="center">Identifiability of word models</h2>
Most biological signaling pathways and reaction systems are present as raw text and in articles which are a common form of communication among the research community. Unfortunately, the natural language descriptions are not computationally viable. Hence this project is an attempt to convert the natural language descriptions into executable models. Furthermore, we also do an identifiability analysis in DAISY from the output generated by the assembled model.



<h3> INDRA model assembly </h3>
Assembling of the model into DAISY amenable format directly from the publication ID or raw text
<h4> Usage  </h4>
<b>Note: You need to run this in Powershell </b>
1. Clone the git repository and run the setup file to install all the dependencies.
```
git clone https://github.com/dhruvi2k/Dual-Degree-Project.git
./ez.setup.ps1
```

2. Change the directory to INDRA_model
```
cd INDRA_model
```

3. To run the INDRA assembly code, use the following command,in case the policy type is passed as a dictionary i.e. different and specific policies for different events
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy '{\"Phosphorylation\":\"two_step\",\"Dephosphorylation\":\"two_step\"}'
```
Use the followint command, in case the policy type is global and has to be passed as a single string. For example "one_step", "two_step", "michaelis_menten"
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy "one_step"
```


4. Two output files will be saved in the output folder. output.txt specifies the desired dynamical model in the DAISY format while attributes.txt mentions the attributes of the assembled model like statements, rules, monomers and the number of ODEs in the reaction system. In case you want to change the name of the output files, you can specify them as follows
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy "one_step" --output_txt "new_output_name.txt" --attr_txt "new_attribute_name.txt"
```

<h4> Tested examples </h4>

1. Simple sentence - MAP21 that is phosphorylated on S218 and S222 phosphorylates MAPK1 at T185. Different policies can be tried
```
python main.py --id_types "raw text" --ids "MAP21 that is phosphorylated on S218 and S222 phosphorylates MAPK1 at T185" --policy "one_step"
```
2. Set of sentences without context - RAF binds Vemurafenib, RAF phosphorylates MEK. MEK phosphorylates ERK
```
python main.py --id_types "raw text" --ids "RAF binds Vemurafenib, RAF phosphorylates MEK. MEK phosphorylates ERK" --policy "one_step"
```
3. Set of sentences with context - RAF binds Vemurafenib. RAF not bound to vemurafenib phosphorylates MEK, Phosphorylated MEK not bound to RAF phosphorylates ERK
```
python main.py --id_types "raw text" --ids "RAF binds Vemurafenib. RAF not bound to vemurafenib phosphorylates MEK, Phosphorylated MEK not bound to RAF phosphorylates ERK" --policy "one_step"
```

4. Given PMCID - 4916225
```
python main.py --id_types "PMCID" --ids "4916225" --policy "one_step"
```

5. Set of publications PMCID - 2132449 PMID - 20668238 PMID 23153539
```
python main.py --id_types "PMCID" "PMID" "PMID" --ids "2132449" "20668238" "23153539" --policy '{\"Phosphorylation\":\"two_step\",\"Dephosphorylation\":\"two_step\"}'
```


<h4> Sending the output file to DAISY and REDUCE </h4>
Note: This will require installation of [DAISY](https://daisy.dei.unipd.it/)

1. Copy the path of the output.txt file generated in the first routine (say `C:\MOD\output.txt`)

2. Run the following command in REDUCE
```
load daisy$
```

3. In order to save results in a file, e.g. daisy_results.txt with the path `C:\MOD\daisy_results.txt`  run the command

```
OUT "C:\MOD\daisy_results.txt"$ 
IN "C:\MOD\output.txt"$ 
SHUT "C:\MOD\daisy_results.txt"$

```

4. Or if you wish to see the results on the REDUCE terminal itself, simply run the command:

```
IN "C:\MOD\output.txt"$ 
```

<h3> Running REACH </h3>

Tags the named entities and identifies biological events using rule based extraction.

<b>Note: You can either run this on a Linux system or Google Colab. </b>

1. Clone the git repository and cd into it
```
git clone https://github.com/clulab/reach.git
cd reach
```

2. Run the following command to install sbt
```
curl -fL https://github.com/coursier/launchers/raw/master/cs-x86_64-pc-linux.gz | gzip -d > cs && chmod +x cs && ./cs setup
```
3. Add the directory `~/.local/share/coursier/bin`  to PATH manually or run this command if in Google Colab
```
import os
os.environ["PATH"] = os.environ["PATH"] + ":~/.local/share/coursier/bin"
print(os.environ["PATH"])
```

4. Build the code
```
sbt compile
```

5. Create the following directory and cd into it
```
mkdir ~/Documents/reach/papers
cd ~/Documents/reach/papers
```
6. Download the fetch python file which will fetch the publication as an nxml file and fetch the publication with a given PMCID as shown below. After the nxml file has been fetched, remove the fetch python file
```
wget https://gist.githubusercontent.com/myedibleenso/f233359445461a71ad37017393fe921f/raw/982275ad8d5070e8c0bc5c07edcfec1cd804c611/fetch_nxml.py
python3 fetch_nxml.py --pmcids PMCID1 PMCID2
rm fetch_nxml.py
```

7. Change the directory to back to the reach folder
```
cd reach
```

8. Run the .sh file will run REACH which will create separate directories for the publications 
```
./runReachCLI.sh 
```
9. Outputs will be stored inside directory `~/Documents/reach/output` with different outputs for fetched publications under the corresponding PMCID. Each directory has different json and text files denoting the extracted named entities and relations.


