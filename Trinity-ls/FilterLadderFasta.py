
# coding: utf-8

# In[1]:


from Bio import SeqIO
import sys
import numpy as np
import bisect


# In[132]:


#inputFastaFile = '/algbio/shounak/raw/simulated/human/newBetas/ImperfectSep100High/Assemblies/denovo/IndividualBand/LadderMerge/trinity_4/trinity_filtered.fasta'
inputFastaFile = sys.argv[1]

#outputFastaFile = '/algbio1/shounak/raw/simulated/human/newBetas/ImperfectSep100High/Assemblies/denovo/LadderMerge/ladderMerge/final.fasta'
outputFastaFile = sys.argv[2]

# filterType : length or id
#filterType = 'length'
filterType = sys.argv[3]



if(filterType == "length"):
    # name of the band where the current fasta file was assembled
    #bandName = "band3"
    bandName = sys.argv[4]
    # Extent of filter
    # EQUAL, CORRECT_PLUS_ONE
    # Required only for Low res filter
    #filterExtent = "EQUAL"
    filterExtent = sys.argv[5]
    # Filter resolution
    # LOW or HIGH
    #filterResolution = "LOW"
    filterResolution = sys.argv[6]
    # Band probability file
    # high res or low res band prob file (only tab separated files)
    #bandProb = np.loadtxt(open("/gcm-lfs1/shounak/LadderSeq/Quantification/exponentialIncreaseBetas/estimatedBetas/simulations/sim_30000000/sim_1/betas_realistic/estimatedBeta_low_text.txt", "rb"), delimiter="\t")
    #bandProb = np.loadtxt(open(sys.argv[7], "rb"), delimiter="\t")
    bandProb = np.genfromtxt(sys.argv[7], delimiter="\t", dtype='str', usecols=(0,1,2,3,4,5,6))
    bandProb = bandProb.astype(np.float)

elif(filterType == "id"):
    #transcriptList = '/algbio1/shounak/raw/simulated/human/newBetas/ImperfectSep100High/Assemblies/denovo/LadderMerge/ladderMerge/transcripts_tobeKept.txt'
    transcriptList = sys.argv[4]

# In[133]:


def getAcceptedBands(transcriptLength,bandList):
    subBandInterval = bisect.bisect_right(bandList,transcriptLength)
    #print(subBandInterval)
    #print(bandProb[subBandInterval-1])

    maxProb = (np.sort(bandProb[subBandInterval-1]))[len(bandProb[subBandInterval-1])-1]
    secondMaxProb = (np.sort(bandProb[subBandInterval-1]))[len(bandProb[subBandInterval-1])-2]
    secondMax_max_prob_ratio = secondMaxProb/maxProb

    if(secondMax_max_prob_ratio>= 0.5):
        acceptedBands = [np.where(bandProb[subBandInterval-1] == maxProb)[0][0],np.where(bandProb[subBandInterval-1] == secondMaxProb)[0][0]]
    else:
        acceptedBands = [np.where(bandProb[subBandInterval-1] == maxProb)[0][0]]

    return acceptedBands



# In[3]:


fasta_sequences = SeqIO.parse(open(inputFastaFile),'fasta')

bandListHigh = (0,910,1083,1268,1440,1590,1757,1920,2070,2273,2462,2701,2932,3230,3554,3893,4525,5324,6378,100000000000)
bandListLow = (0,1000,1500,2000,3000,4000,6000,100000000000)

if(filterType == 'length'):
    with open(outputFastaFile, "w") as handle:
        for record in fasta_sequences:
            transcriptLength = len(record)
            if(filterResolution == "HIGH"):
                acceptedBands = getAcceptedBands(transcriptLength,bandListHigh)
            elif(filterResolution == "LOW"):
                acceptedBands = getAcceptedBands(transcriptLength,bandListLow)

            if(filterExtent == "CORRECT_PLUS_ONE"):
                acceptedBands.append(max(acceptedBands)+1)

            if("band" in bandName):
                # subtracting 1 here since np arrays have a 0 based index whereas bands have a one based index
                assembledBand = int(bandName[-1:])-1
                if(assembledBand in acceptedBands):
                    SeqIO.write(record, handle, "fasta")
            elif("comb" in bandName):
                # subtracting 1 here since np arrays have a 0 based index whereas bands have a one based index
                assembledBand1 = (int(bandName[-1:]) -1) - 1
                assembledBand2 = (assembledBand1 + 1)
                if((assembledBand1 in acceptedBands) or (assembledBand2 in acceptedBands)):
                    SeqIO.write(record, handle, "fasta")

elif(filterType == 'id'):
    transcriptList = [line.rstrip('\n') for line in open(transcriptList)]
    with open(outputFastaFile, "w") as handle:
        for record in fasta_sequences:
            transcriptId = record.id
            if(transcriptId in transcriptList):
                SeqIO.write(record, handle, "fasta")
