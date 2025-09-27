# Mouse

## B6

### mm39

### B6 T2T

UCSC Genome browser: [**C57BL_6J_T2T_v1 Jun. 2024 Genome Browser** - GCA_964188535.1 assembly](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=3166850614_KGl5M3PHqx7rQus5Bf8PDXS9WRjS)

#### Repeats

Repeat Masker annotations on the GCA_964188535.1_C57BL_6J_T2T_v1 genome assembly from [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?db=hub_6517761_GCA_964188535.1&hgta_group=varRep&hgta_track=hub_6517761_repeatMasker&hgta_table=hub_6517761_repeatMasker&hgta_doSchema=describe+table+schema)
Get chromosome names from:   [GCA_964188535.1.chromAlias.txt](https://hgdownload.soe.ucsc.edu/hubs/GCA/964/188/535/GCA_964188535.1/GCA_964188535.1.chromAlias.txt) from this [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=3166850614_KGl5M3PHqx7rQus5Bf8PDXS9WRjS)

```bash
# download repeats
wget https://hgdownload.soe.ucsc.edu/gbdb/genark/GCA/964/188/535/GCA_964188535.1/bbi/GCA_964188535.1_C57BL_6J_T2T_v1.rmsk.bb
# covert from bigBed to bed
bigBedToBed GCA_964188535.1_C57BL_6J_T2T_v1.rmsk.bb GCA_964188535.1_C57BL_6J_T2T_v1.rmsk.bed

# download chromosome names
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/964/188/535/GCA_964188535.1/GCA_964188535.1.chromAlias.txt

# change chromosome names with script
bash map_chromosomes_for_b6_repeats.sh GCA_964188535.1.chromAlias.txt GCA_964188535.1_C57BL_6J_T2T_v1.rmsk.bed C57BL_6J_T2T_v1_repeats.bed
```

```bash fold title:"map_chromosomes_for_b6_repeats.sh"
chrom_names=$1
input_file=$2
output_file=$3

awk '
    # This block runs ONLY for the chromosome names file ($chrom_names)
    BEGIN { OFS="\t" }
    FNR==NR {
        # RULE 1: If the identifier starts with "CAX"
        if ($1 ~ /^CAX/) {
            # Use the entire second column as the new name
            map[$1] = $2
        }
        # RULE 2: If the identifier starts with "OZ"
        else if ($1 ~ /^OZ/) {
            # Split the second column by "#" and get the last part
            n = split($2, parts, "#")
            map[$1] = parts[n]
        }
        next
    }

    # This block runs for the repeat file ($input_file)
    {
        # Check if the current line ID is in our map
        if ($1 in map) {
            # If it is, replace it with the new name
            $1 = map[$1]
        }
        # Print the line (modified or original)
        print
    }
' $chrom_names $input_file > $output_file
```


## CAST 

### CAST T2T

UCSC Genome browser: [**CAST_EiJ_T2T_v1 Jun. 2024 Genome Browser** - GCA_964188545.1 assembly](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=3166850614_KGl5M3PHqx7rQus5Bf8PDXS9WRjS)

#### Repeats

Repeat Masker annotations on the GCA_964188545.1_CAST_EiJ_T2T_v1 genome assembly from [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?db=hub_6552617_GCA_964188545.1&hgta_group=varRep&hgta_track=hub_6552617_repeatMasker&hgta_table=hub_6552617_repeatMasker&hgta_doSchema=describe+table+schema)
Get chromosome names from:   [GCA_964188545.1.chromAlias.txt](https://hgdownload.soe.ucsc.edu/hubs/GCA/964/188/545/GCA_964188545.1/GCA_964188545.1.chromAlias.txt) from this [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=3166850614_KGl5M3PHqx7rQus5Bf8PDXS9WRjS)

```bash
# download repeats
wget https://hgdownload.soe.ucsc.edu/gbdb/genark/GCA/964/188/545/GCA_964188545.1/bbi/GCA_964188545.1_CAST_EiJ_T2T_v1.rmsk.bb
# covert from bigBed to bed
bigBedToBed GCA_964188545.1_CAST_EiJ_T2T_v1.rmsk.bb GCA_964188545.1_CAST_EiJ_T2T_v1.rmsk.bed

# download chromosome names
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/964/188/545/GCA_964188545.1/GCA_964188545.1.chromAlias.txt

# change chromosome names with script
bash map_chromosomes_for_cast_repeats.sh GCA_964188545.1.chromAlias.txt GCA_964188545.1_CAST_EiJ_T2T_v1.rmsk.bed CAST_EiJ_T2T_v1_repeats.bed
```

```bash fold title:"map_chromosomes_for_cast_repeats.sh"
chrom_names=$1
input_file=$2
output_file=$3

awk '
    # This block runs ONLY for the chromosome names file ($chrom_names)
    BEGIN { OFS="\t" }
    FNR==NR {
        # RULE 1: If the identifier starts with "CAX"
        if ($1 ~ /^CAX/) {
            # Split the second column by "." and get the last part
            n = split($2, parts, /\./)
            map[$1] = parts[n]
        }
        # RULE 2: If the identifier starts with "OZ"
        else if ($1 ~ /^OZ/) {
            # Split the second column by "#" and get the last part
            n = split($2, parts, "#")
            map[$1] = parts[n]
        }
        next
    }

    # This block runs for the repeat file ($input_file)
    {
        # Check if the current line ID is in our map
        if ($1 in map) {
            # If it is, replace it with the new name
            $1 = map[$1]
        }
        # Print the line (modified or original)
        print
    }
' $chrom_names $input_file > $output_file
```