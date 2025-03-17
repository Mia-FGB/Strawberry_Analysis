#!/usr/bin/python

import sys, getopt, errno

#A class which contains all the info for eqch alignment 
class alignment_info:
    def __init__(self, q_name, taxaID, match_bases, map_length, MQ):
        self.q_name = q_name
        self.taxaID = taxaID
        self.match_bases = match_bases
        self.map_length = map_length
        self.MQ = MQ
    #Function within class that will print query name & MQ 
    def print_alignment(self):
        print(self.q_name, self.MQ)

#This block of code means it can take any paf file 
# When using the -i flag in the terminal use the filename given as input_arg
try:
    opts, args = getopt.getopt(sys.argv[1:],"i:")
except getopt.GetoptError:
    print("Option not recognised.")
    print("python my_script.py -i <input_argument>")
    sys.exit(2)
for opt, arg in opts:
    if opt == "-i":
        input_arg = arg
    else:
        print("python my_script.py -i <input_argument>")
        sys.exit(2)

#From here on can use pafFilename rather than input_arg
pafFilename = input_arg
#setting up a dictionary with queries as keys 
queries = dict()

#A dictionary which wil have taxaID for keys and counts as values 
taxa_count = dict()

#parsing in the paf file 
try:
    with open(pafFilename, 'r') as pafFile: 
        for line in pafFile:
            line = line.strip()
            fields = line.split()
            q_name = fields[0].strip()
            #only want to get the numbers not the whole section
            taxaID =  fields[5].split("|")[1]
            match_bases = int(fields[9].strip())
            map_length = int(fields[10].strip())
            MQ = int(fields[11].strip())
            #starting the dictionary to populate with taxaIDs
            taxa_count[taxaID] = 0
            #only want to include MQ which are 0 or >=5 so filter on this line before creating the dictionary 
            if MQ >= 5 or MQ == 0:
                #Now calling the class alignment_info and inputting this info
                ai = alignment_info(q_name, taxaID, match_bases, map_length, MQ)
                #whilst looping through each line in the paf file
                if q_name not in queries.keys():
                    #adding a new key for new queries
                    queries[q_name] = list()
                #then appending the class to the dictionary based on queryName
                queries[q_name].append(ai)
except (OSError, IOError) as e: 
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print ("Could not find file " + pafFilename)
		sys.exit(2)

#Function to select which alignments are retained
def parse(queries):
    single_query = []
    multiple_query = []
    #counting how many reads are ignored
    ignored_reads = 0
    #looping through the queries dictionary 
    for q_name in queries.keys():
        #creating a list that contains the info from the dictionaty 
        ai_list = queries[q_name]
        #If there is only one input for a query name then the read is unique 
        if len(ai_list) == 1:
            #so it can be added to the unique list 
            single_query.append(q_name)
            #this ai list only contains one thing so we can get it with [0]
            ai = ai_list[0]
            #getting the taxaID out of the class we set up using the ai_list query name
            taxa_id = ai.taxaID
            #updating the count for the specific taxaID in the dictionary
            #taxa_ID is the key
            taxa_count[taxa_id] += 1
        else:
            #If there is more than one entry then there must be multiple reads
            multiple_query.append(q_name)
            #making a set as this can only contain unique values 
            mapping_qualities = set()
            multiple_taxa_ID = set()
            #looping through list of dict info
            for ai in ai_list:
                #adding the mapping quality to the set
                mapping_qualities.add(ai.MQ)
                multiple_taxa_ID.add(ai.taxaID)

            #1x taxa, add as 1 alignment
            if len(multiple_taxa_ID) == 1:
                #extracting from the set to use it in dict
                taxa_id = list(multiple_taxa_ID)[0]
                taxa_count[taxa_id] +=1
            
            #1xMQ multi x taxa
            #currenlty going to ignore these reads but want to print them to see how big of a problem this is 
            #In the future could look at LCA 
            if len(mapping_qualities) == 1 and len(multiple_taxa_ID) > 1:
                #print("This query: " , ai.q_name , "has the same MQ", mapping_qualities, "but maps to different taxaIDs: ", multiple_taxa_ID, " so has been ignored.")
                ignored_reads += 1

            #multi x MQ and multi x taxaID, just want to take taxaID from highest MQ
            if len(mapping_qualities) > 1 and len(multiple_taxa_ID) > 1:
                #sort and choose the highest MQ to add to taxa_count dictionary
                maxMQ = 0
                taxaIDForMaxMQ = ""
                for ai in ai_list:
                    #looping though the ai list which contains all info
                    #Still within the q_name loop
                    if ai.MQ > maxMQ:
                        maxMQ = ai.MQ
                        taxaIDForMaxMQ = ai.taxaID
                        taxa_count[taxaIDForMaxMQ] +=1
   
    #prints the dictionary as a list of taxaIDs & counts
    for taxa_id, count in taxa_count.items():
        print('{}' '\t' ' {}'.format(taxa_id, count))

    print('This many reads are ignored: {} {} '.format(ignored_reads, input_arg))

    




print(parse(queries))


