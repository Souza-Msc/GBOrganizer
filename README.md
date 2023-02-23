# GB Organizer

Genbank is the biggest database of genetic/genomic material in the world. The website is used every day by a large part of the scientific community, either for depositing data or for collecting it. They receive all types of genomic and genetic data from all living organisms from every part of the globe. Although useful, this huge flow of information comes with a price. They don't have a high level of filtering of the data entered into their database, so the first step for all scientists is to sort out the material they intend to use. This process is usually massive. To help the students in my lab with this first step, I created this script. It'll sort the genetic/genomic data and compile the important information of each sequence in a dataframe, which is much easier to use.

To utilize this script, all you need to do is download the .gb file from the Genbank website (https://www.ncbi.nlm.nih.gov/genbank/) and have Python installed on your device. You can use this script to compile only the information you need in a dataframe by using the following command:

```
python GBOrganizer.py --input_file sequence.gb
```

If you want the code to split the data into fasta files containing the same type of marker, use:

```
python GBOrganizer.py --input_file sequence.gb --split_file Yes
```

You can also set a minimum length to further specify your sorting. To do this, use:

```
python GBOrganizer.py --input_file sequence.gb --split_file Yes --lenght 1000

## Contact me

Questions, suggestions, comments, etc? Just send a e-mail to pedromsouza0@gmail.com or souza.pedro@ecologia.ufjf.br
