import requests, sys, json, pickle, nltk, random, operator, traceback, os
import xml.etree.ElementTree as etree
from nltk.corpus import stopwords
from misc_functions import loadJson, writeJson, printList, fileToList, listToFile, getKey
from time import sleep
from collections import OrderedDict

# Input Files
COMPOUND_DICTIONARY_FILE = "Input\\compound_dictionary.json"
SMILES_DICTIONARY_FILE = "Input\\smiles_dictionary.json"
WITHOUT_DETAILS_FILE = "Input\\pdbs_without_details.json"
STOP_WORDS_FILE = "Input\\stop_words.json"
UNKNOWN_LIST_FILE = "Input\\unknown_list.json"
ERROR_LIST_FILE = "Input\\error_list.json"
LOWERCASE_REPLACEMENT_FILE = "Input\\replacementLowercase.json"
SENSITIVE_REPLACEMENT_FILE = "Input\\replacementSensitive.json"
MIXTURES_FILE = "Input\\mixture_compounds.json"

STRUCTURES_FILE = "Structures\\structures.pkl"
SENSIBLE_STRUCTURES_FILE = "Structures\\sensible_structures.pkl"

if not os.path.exists("Output"):
    os.makedirs("Output")

# Output Files
DETAILS_FILE = "Output\\details.txt"
SENSIBLE_DETAILS_FILE = "Output\\sensible_details.txt"
ERROR_DETAILS_FILE = "Output\\error_details.txt"
COMPOUND_FREQUENCY_FILE = "Output\\compound_frequency.txt"
UNKNOWN_FREQUENCY_FILE = "Output\\unknown_frequency.txt"
PENDING_FREQUENCY_FILE = "Output\\pending_frequency.txt"

# Load lists and dictionaries from files
print("Loading dictionary files...")
try:
    unknownList = loadJson(UNKNOWN_LIST_FILE)
    compoundDictionary = loadJson(COMPOUND_DICTIONARY_FILE)
    smilesDictionary = loadJson(SMILES_DICTIONARY_FILE)
    errorList = loadJson(ERROR_LIST_FILE)
    lowercaseReplacement = loadJson(LOWERCASE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    sensitiveReplacement = loadJson(SENSITIVE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    mixturesDictionary = loadJson(MIXTURES_FILE)
    EXTRA_STOP_WORDS = set(loadJson(STOP_WORDS_FILE))
    STOP_WORDS = set(stopwords.words("english")) - {'m', 'am'} | EXTRA_STOP_WORDS - {'m'}
except FileNotFoundError as notFound:
    print(notFound)
    print("ERROR: The file {} cannot be found. Verify that it is in the proper directory.".format(notFound.filename))

# Compounds that contain numbers (eg jeffamine 600)
NUMBERED_COMPOUNDS = ["jeffamine", "propoxylate"]

class Structure:

    # Global variables used for parsing of details
    RESEVOIR_INDICATOR_WORDS = {"reservoir", "resevoir", "crystallization", "precipitant", "well", "solution", "cocktail"}


    def __init__(self, pdbid, pmcid, details, compounds, pH, temperature, method, sequences, resolution):
        self.pdbid = pdbid # String
        self.pmcid = pmcid # String
        self.details = details # crystallization details, String
        self.compounds = compounds # List of compounds followed by concentration
        self.pH = pH # float
        self.temperature = temperature # float, in K
        self.method = method # crystallization method, String
        self.sequences = sequences # List
        self.resolution = resolution # float, in angstroms
        self.compounds = []
        #self.parseDetails() # Sets the compounds list

    def __str__(self):
        return("""
        PDB: {0}
        PMCID: {1}
        Compounds: {2}
        pH: {3}
        Temperature: {4}
        Method: {5}
        Resolution: {6}
        """.format(self.pdbid, self.pmcid, self.compounds, self.pH, self.temperature, self.method, self.resolution))

    def parseDetails(self, debug=False):
        """Parses the details string and returns a list of compounds, followed by concentration, or None if conc. is not found
        Format of compounds = ['compound name', '100', 'another compound name', '45%']
        If debug=True, then this function will print the words after every step
        Also sets the compounds list
        Returns None if failed or if details are unavailible
        """

        if self.details == None:
            return None

        details = self.details

        if debug:
            print("Debugging PDB {}...\n".format(self.pdbid))
            print("Raw details:\n"+details+"\n")

        # Add spaces after commas before numbers (eg "sodium acetate,4% PEG4K" --> "sodium acetate, 4% PEG4K")
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(details)-1):
                if details[j] == "," and details[j-1].isalpha() and isNumber(details[j+1]):
                    # Insert a space
                    details = details[:j+1] + " " + details[j+1:]
                    runAgain = True
                    break

        if debug:
            print("Dealing with lack of space after commas:\n"+details+"\n")

        details = wordReplacement(details, sensitiveReplacement, debug=debug) # Replace with case sensitivity
        details = wordReplacement(details.lower(), lowercaseReplacement, debug=debug) # Replace without case sensitivity

        # Remove all text after 'cryo'
        cryoLocation = details.find("cryo")
        if cryoLocation != -1:
            details = details[:cryoLocation]

        if debug:
            print("Word replacement:\n"+details+"\n")

        words = nltk.word_tokenize(details)

        if debug:
            print("Tokenized words:\n"+str(words)+"\n")

        # Get correct sentence
        if len(words) > 1 and words[0] == "protein" and words[1] in Structure.RESEVOIR_INDICATOR_WORDS:
            reservoirSolutionFound = False
            for j in range(0, len(words)-1):
                if ((words[j] in Structure.RESEVOIR_INDICATOR_WORDS) and (words[j+1] in Structure.RESEVOIR_INDICATOR_WORDS)):
                    words = words[j:]
                    reservoirSolutionFound = True
                    break
            if (not reservoirSolutionFound):
                words = ["ERROR", "Resevoir solution not found"]

        if debug:
            print("Correct sentence:\n"+str(words)+"\n")

        # Fix commas between numbers (turn 6000,10 in to 6000 , 10)
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if "," in words[j] and isNumber(words[j].replace("%","").replace("m","").replace("k","")):
                    numbers = words[j].split(",")
                    if len(numbers) == 2 and (len(numbers[1]) != 3 or len(numbers[0]) > 3 or "." in numbers[1]) and numbers[0] != "" and numbers[1] != "":
                        runAgain = True
                        words[j] = numbers[0]
                        words.insert(j+1, ",")
                        words.insert(j+2, numbers[1])
                        if len(numbers[1]) > 2 and numbers[1][-2:] == "mm":
                            words[j+2] = numbers[1][:-2]
                            words.insert(j+3, "mm")
                        elif numbers[1][-1:] == "m" or numbers[1][-1:] == "%":
                            words[j+2] = numbers[1][:-1]
                            words.insert(j+3, numbers[1][-1])
                        break

        if debug:
            print("Fix commas:\n"+str(words)+"\n")

        # Process micromolar (uM) concentrations, add to number like: "35 um" --> "35um"
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(1, len(words)):
                if words[j] == "um" and isNumber(words[j-1]):
                    words[j-1] = words[j-1]+"um"
                    del words[j]
                    runAgain = True
                    break

        if debug:
            print("Process uM concentrations:\n"+str(words)+"\n")

        # Remove pH and temperature
        runAgain = True
        while (runAgain): # pH
            runAgain = False
            for j in range(0, len(words)):
                if j+1 < len(words) and words[j] == "ph" and isNumber(words[j+1]):
                    del words[j:j+2]
                    runAgain = True
                    break
                elif len(words[j]) > 2 and words[j][:2] == "ph" and isNumber(words[j][2:].replace("=","").replace(" ","")):
                    del words[j]
                    runAgain = True
                    break
                elif j+2 < len(words) and words[j] == "ph" and words[j+1] == "=" and isNumber(words[j+2]):
                    del words[j:j+3]
                    runAgain = True
                    break
        runAgain = True
        while (runAgain): # Temperature
            runAgain = False
            for j in range(1, len(words)):
                if (words[j][-1] == 'k' and len(words[j]) > 3 and isNumber(words[j][:-1])):
                    del words[j]
                    runAgain = True
                    break
                elif words[j] == 'k' and isNumber(words[j-1]):
                    temp = float(words[j-1].replace(",",""))
                    if temp > 200 and temp < 400:
                        del words[j-1:j+1]
                        runAgain = True
                        break

        if debug:
            print("Remove pH and temperature:\n"+str(words)+"\n")

        # Process compounds with numbers in them
        for j in range(0, len(words)-1):
            if words[j] in NUMBERED_COMPOUNDS and isNumber(words[j+1]):
                words[j] = words[j] + " " + words[j+1]
                del words[j+1]
                break

        # Process PEG and MPEG compounds
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if (words[j][:3] == "peg" or words[j][:4] == "mpeg"):
                    pegNumber = ""
                    if j+1<len(words) and words[j+1] in STOP_WORDS:
                        runAgain = True
                        del words[j+1]
                        break
                    elif len(words[j]) > 3 and words[j][:3] == "peg": # If the number is already attached, eg PEG5000
                        if j+1<len(words) and words[j+1] == "mme": # Ex. PEG3000 MME
                            pegNumber = "MME " + words[j][3:].replace(",","").replace("-", " ")
                            del words[j+1]
                        else:
                            pegNumber = words[j][3:].replace("-", " ").replace(",","")
                    elif len(words[j]) > 4 and words[j][:4] == "mpeg": #MPEG4000
                        pegNumber = words[j][4:].replace("-", " ").replace(",","")
                    elif j+2 < len(words) and (words[j+1] == "mme" or words[j+1] == "monomethyl") and words[j+2] not in STOP_WORDS: # ex. PEG MME 5000
                        if isNumber(words[j+2].replace("k","000")):
                            pegNumber = "MME " + words[j+2].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+2 < len(words) and words[j+2] == "mme" and words[j+1] not in STOP_WORDS: # ex. PEG 5000 MME
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = "MME " + words[j+1].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+1<len(words) and words[j+1] not in STOP_WORDS: # Ex PEG 5000
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = words[j+1].replace(",","")
                            del words[j+1]
                    if len(pegNumber) > 3 and pegNumber[-3:] == "mme":
                        pegNumber = "MME " + pegNumber[:-3]
                    pegNumber = pegNumber.replace("000k","000").replace("k","000").replace("-","")
                    if words[j][:4] == "mpeg":
                        words[j] = "MPEG "+pegNumber
                    else:
                        words[j] = "PEG "+pegNumber
                    runAgain = True
                    break

        if debug:
            print("Process PEG and MPEG compounds:\n"+str(words)+"\n")

        # Remove terminal periods from numbers
        for j in range(0, len(words)):
            if isNumber(words[j]) and words[j][-1] == ".":
                words[j] = words[j][:-1]

        if debug:
            print("Remove terminal periods from numbers:\n"+str(words)+"\n")

        # Put together numbers separated by "-"
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(words)-1):
                if words[j] == "-" and isNumber(words[j-1]):
                    if isNumber(words[j+1]) or (words[j+1][-1] == "m" and isNumber(words[j+1][:-1])) or (len(words[j+1]) > 2 and words[j+1][-2] == "mm" and isNumber(words[j+1][:-2])):
                        words[j-1] = printList(words[j-1:j+2], "")
                        del words[j:j+2]
                        runAgain = True
                        break

        if debug:
            print("Combine numbers separated by '-':\n"+str(words)+"\n")

        # Average ranges of numbers
        for j in range(0, len(words)):
            if '-' in words[j]:
                if averageNumberString(words[j]) != None: # If there is no letter in front of the numbers
                    words[j] = str(averageNumberString(words[j]))
                elif averageNumberString(words[j][:-1]) != None: # If there is one letter in front of the numbers (e.g. 'm' or '%')
                    words.insert(j+1, words[j][-1])
                    words[j] = str(averageNumberString(words[j][:-1]))
                elif averageNumberString(words[j][:-2]) != None: # If there are two letters in front of the numbers (i.e. Mm)
                    words.insert(j+1, words[j][-2:])
                    words[j] = str(averageNumberString(words[j][:-2]))

        if debug:
            print("Average ranges of numbers:\n"+str(words)+"\n")

        if len(words) == 0:
            return None

        # Process concentrations
        runAgain = True
        while(runAgain):
            for j in range(0, len(words)):
                if words[j][-1] == "m" and isNumber(words[j][:-1]): # Units are M
                    words[j] = str(float(words[j][:-1].replace("%","").replace(",",""))*1000)
                    runAgain = True
                    break
                elif words[j][-2:] == "mm" and isNumber(words[j][:-2]): # Units are mM
                    words[j] = words[j][:-2]
                    runAgain = True
                    break
                elif words[j][-2:] == "um" and isNumber(words[j][:-2]): # Units are uM
                    words[j] = str(float(words[j][:-2].replace("%","").replace(",",""))/1000)
                    runAgain = True
                    break
                elif isNumber(words[j]) and j+1 < len(words): # Conc. is separated from number by space, e.g. "0.1 M"
                    if words[j+1] == "m":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",",""))*1000)
                        runAgain = True
                        break
                    elif words[j+1] == "mm":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",","")))
                        runAgain = True
                        break
                    elif words[j+1] == "%":
                        if j+2 < len(words) and (words[j+2] == "w/v" or words[j+2] == "v/v"):
                            units = "% " + words[j+2]
                            del words[j+2]
                        else:
                            units = "%"
                        del words[j+1]
                        words[j] = str(words[j]+units)
                        runAgain = True
                        break
                else:
                    runAgain = False

        if debug:
            print("Process concentrations:\n"+str(words)+"\n")

        # Remove excess punctuation
        for j in range(0, len(words)):
            if words[j] in {",", "."}:
                if all(word in {",", "."} for word in words[j:]):
                    del words[j:]
                    break

        if debug:
            print("Remove excess punctuation:\n"+str(words)+"\n")

        # Combline parts of chemical compounds
        runAgain = True
        while runAgain:
            runAgain = False
            for j in range(0, len(words)-1):
                if isCompound(words[j]) and isCompound(words[j+1]) and words[j][:4] != "PEG " and words[j+1][:4] != "PEG " and words[j][:4] != "MPEG " and words[j+1][:4] != "MPEG ":
                    k = j
                    chemicalCompound = []
                    while(k < len(words) and isCompound(words[k])):
                        chemicalCompound.append(words[k])
                        k += 1
                    del words[j+1:k]
                    words[j] = printList(chemicalCompound, " ")
                    runAgain = True
                    break

        if debug:
            print("Combine parts of chemical compounds:\n"+str(words)+"\n")

        # Remove periods from compounds
        for j in range(0, len(words)):
            if isCompound(words[j]):
                words[j] = words[j].replace(".","")

        if debug:
            print("Remove periods from compounds:\n"+str(words)+"\n")

        # Combine w/v or v/v with non-percent nubers, and assume it's supposed to be percent
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)-1):
                if isNumber(words[j]) and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                    words[j] = words[j] + "% " + words[j+1]
                    del words[j+1]
                    runAgain = True
                    break

        if debug:
            print("Combine w/v and v/v with non-percent numbers:\n"+str(words)+"\n")

        # Determine if concentration or compounds come first
        concentrationBeforeCompound = True
        for j in range(0, len(words)-1):
            if isConcentraton(words[j]) and isCompound(words[j+1]):
                concentrationBeforeCompound = True
                break
            elif isConcentraton(words[j+1]) and isCompound(words[j]):
                concentrationBeforeCompound = False
                break

        # Remove stop words
        words = [w for w in words if w not in STOP_WORDS]

        if debug:
            print("Remove stop words:\n"+str(words)+"\n")

        # Add compounds to list
        compounds = []
        if len(words) > 0 and words[0] != "ERROR":
            runAgain = True
            while(runAgain):
                runAgain = False
                for j in range(0, len(words)):
                    if isCompound(words[j]):
                        if concentrationBeforeCompound and j > 0 and isConcentraton(words[j-1]): # If a concentration is found (before compound)
                            compound = words[j]
                            concentration = words[j-1]
                            # If w/v or v/v comes after compound
                            if j+1<len(words) and "%" in concentration and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                                concentration = concentration + " " + words[j+1]
                                del words[j-1:j+2]
                            else:
                                del words[j-1:j+1]
                        elif not concentrationBeforeCompound and j+1 < len(words) and isConcentraton(words[j+1]): # if a concentration is found (after compound)
                            compound = words[j]
                            concentration = words[j+1]
                            del words[j:j+2]
                        else: # No concentration found
                            compound = words[j]
                            concentration = None
                            del words[j]
                        compounds.append(compound)
                        compounds.append(concentration)
                        runAgain = True
                        break
        self.compounds = compounds

        if debug:
            print("Extract compounds:\n"+str(compounds)+"\n")

        return compounds

    def standardizeNames(self):
        """Standardizes the name of an individual structure"""
        runAgain = True
        while(runAgain):
            runAgain = False
            for i in range(0, len(self.compounds), 2):
                if getKey(self.compounds[i]) in compoundDictionary:
                    self.compounds[i] = compoundDictionary[getKey(self.compounds[i])]

                    # Parse multiple compounds delimited with " / "
                    if " / " in self.compounds[i]:
                        multipleCompounds = self.compounds[i].split(" / ")
                        concentration = self.compounds[i+1]
                        for c in multipleCompounds:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                self.compounds.append(concentration)
                            else:
                                print("ERROR: Compound {} was found in multiple compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    self.compounds.append(concentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

                    # Parse mixture compounds, that is,
                    # names of mixtures which refer to multiple chemical compounds
                    if self.compounds[i] in mixturesDictionary:
                        mixtureDict = mixturesDictionary[self.compounds[i]] # Dictionary of the compounds in the mixture
                        concentration = self.compounds[i+1]
                        if concentration != None:
                            findPercent = concentration.find("%")
                            percentUnits = ""
                            if findPercent != -1:
                                percentUnits = concentration[findPercent:] # Everything after and including '%'
                                concentration = concentration[:findPercent] # Everything berfore "%"
                        for c in mixtureDict:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                if concentration != None:
                                    if percentUnits == "":
                                        newConcentration = str(float(concentration)*mixtureDict[c])
                                    else:
                                        # If concentration is a percent, then calculate millimolar concentration
                                        newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                else:
                                    newConcentration = None
                                self.compounds.append(newConcentration)
                            else:
                                print("ERROR: Compound {} was found in mixture compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    if concentration != None:
                                        if percentUnits == "":
                                            newConcentration = str(float(concentration)*mixtureDict[c])
                                        else:
                                            # If concentration is a percent, then calculate millimolar concentration
                                            newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                    else:
                                        newConcentration = None
                                    self.compounds.append(newConcentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

    def hasUnknown(self): # boolean
        """returns True if the structure has a compound in the unknownList"""
        for compound in self.compounds[::2]:
            if compound in unknownList:
                return True

    def hasError(self): # boolean
        """Returns True if the structure has a compound in the errorList"""
        for compound in self.compounds[::2]:
            if compound in errorList:
                return True

    def printError(self, errorMessage, errorObject):
        """Prints an error message regarding a certain structure.
        errorMessage is a string which describes when the error occrued
        errorObject is the exception object which is caught"""
        print("--------------------\nERROR: {} for PDB ID {}.\n".format(errorMessage, self.pdbid))
        traceback.print_tb(errorObject.__traceback__)
        print(str(errorObject)+"\n")
        print("Crystallization Details: {}\n\nCompounds: {}\n--------------------\n".format(self.details, self.compounds))

def parseAllDetails(structureList, searchString=None):
    """Reparses all of the details for a list of structures
    Should be called when the parseDetails function has been modified
    If search string is not None, then it will only parse Structures
    which have the search string in their details, making it faster
    """
    print("Parsing details with search string: \"{}\"...".format(searchString))
    count = 1
    for structure in structureList:
        if count == 1 or count % 10000 == 0:
            print("Parsing structure {} of {}".format(count, len(structureList)))
        try:
            if searchString == None or searchString.lower() in structure.details.lower():
                structure.parseDetails() # void
        except Exception as e:
            structure.printError("Unable to parse details", e)
        count += 1
    print("Writing to file...")
    writeStructures(structureList, STRUCTURES_FILE)

def standardizeAllNames(structureList): # void
    """Standardizes names of a list of compounds based on the compound Dictionary
    Also parses dictionary values which represent multiple compounds (e.g. acetic acid / sodium acetate)
    Also parses mixture compounds, which are mixtures of multiple compounds (e.g Molecular Dimensions Buffers)
    """
    updateDictionary()
    print("Standardizing chemical names...")
    for structure in structureList:
        try:
            structure.standardizeNames()
        except Exception as e:
            structure.printError("Unable to standardize compound names", e)

def updateSmilesDictionary(sensibleStructureList): # void
    """Takes a list of SENSIBLE structures and adds newly
    recognized compounds as blank keys to the SMILES dictionary"""
    updateDictionary()
    print("Updating SMILES dictionary...")
    for structure in sensibleStructureList:
        for compound in structure.compounds[::2]:
            if compound not in smilesDictionary and compound in compoundDictionary.values():
                smilesDictionary[compound] = ""
    writeJson(smilesDictionary, SMILES_DICTIONARY_FILE, indent=2, sort_keys=True)

def getSensibleStructures(structureList):
    """Returns a list of structures that make sense
    That is, all of their compounds are found in the dictionary"""
    standardizeAllNames(structureList)
    print("Searching for sensible structures...")
    sensibleStructures = []
    for structure in structureList:
        isSensible = True
        if structure.compounds == []:
            isSensible = False
        for compound in structure.compounds[::2]:
            if compound not in set(compoundDictionary.values()):
                isSensible = False
                break
        if isSensible:
            sensibleStructures.append(structure)
    return sensibleStructures

def exportSensibleStructures(structureList):
    """Exports a list of sensible structures to a serialized pickle file
    Also exports a details text file with only the sensible structures"""
    sensibleStructures = getSensibleStructures(structureList)
    writeStructures(sensibleStructures, SENSIBLE_STRUCTURES_FILE)
    print("Retrieved {} sensible structures".format(len(sensibleStructures)))

    outputList = []
    print("Exporting structure details...")
    for struc in sensibleStructures:
        outputList.append(struc.details)
        outputList.append(struc.compounds)
    listToFile(outputList, SENSIBLE_DETAILS_FILE)

def exportErrors(structureList, outputFilename=ERROR_DETAILS_FILE, unknown=False):
    """Exports the details of all structures with errors to a text file
    If unknown is True, also export structures with unknown compounds"""
    print("Exporting error details to \"{}\"...".format(outputFilename))
    outputList = []
    for structure in structureList:
        if structure.hasError() or (unknown and structure.hasUnknown()):
            outputList.append(structure.details)
            outputList.append(structure.compounds)
    listToFile(outputList, outputFilename)

def exportAllDetails(structureList, outputFilename=DETAILS_FILE):
    """Exports the details of all structures to a text file"""
    print("Exporting details to \"{}\"...".format(outputFilename))
    outputList = []
    for structure in structureList:
        if structure.compounds != []:
            outputList.append(structure.pdbid + ": " + structure.details)
            outputList.append(structure.compounds)
    listToFile(outputList, outputFilename)

def exportCompoundFrequencies(structureList, outputFilename, mode="recognized"):
    """Takes a list of structures and outputs compound frequency to a txt file
    Mode = "recognized": only export compounds found in the dictionary
    Mode = "unknown": only export compounds found in unknownList
    Mode = "pending": only export compounds in neither dictionary nor unknown list, pending classification
    """
    print("Exporting compound frequencies to \"{}\"...".format(outputFilename))
    outputDictionary = {"Total Compounds": 0}
    if mode == "recognized":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                c = compound
                if getKey(c) in compoundDictionary:
                    c = compoundDictionary[getKey(c)]
                if c in compoundDictionary.values():
                    if c not in outputDictionary:
                        outputDictionary[c] = 1
                    else:
                        outputDictionary[c] += 1
                    outputDictionary["Total Compounds"] += 1

    if mode == "unknown":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                if getKey(compound) in unknownList and getKey(compound) not in compoundDictionary:
                    if compound not in outputDictionary:
                        outputDictionary[compound] = 1
                    else:
                        outputDictionary[compound] += 1
                    outputDictionary["Total Compounds"] += 1

    if mode == "pending":
        for structure in structureList:
            for compound in structure.compounds[::2]:
                if getKey(compound) not in unknownList and getKey(compound) not in compoundDictionary and compound not in compoundDictionary.values():
                    if compound not in outputDictionary:
                        outputDictionary[compound] = 1
                    else:
                        outputDictionary[compound] += 1
                    outputDictionary["Total Compounds"] += 1

    # Export
    outputList = []
    for key, value in sorted(outputDictionary.items(), key=operator.itemgetter(1), reverse=True):
        outputList.append("{:30s}: {:20d}".format(key, value))
    listToFile(outputList, outputFilename)

def updateDictionary():
    """Makes sure all values in the compound dictionary are also keys
    This prevents inconsistency in other functions which assume that this is the case"""
    for value in set(compoundDictionary.values()):
        if getKey(value) not in compoundDictionary:
            compoundDictionary[getKey(value)] = value
    writeJson(compoundDictionary, COMPOUND_DICTIONARY_FILE, indent=2)

def getStructure(structureList, pdbid): # Structure
    """Returns a specific structure in a list based on its pdbid"""
    for structure in structureList:
        if pdbid == structure.pdbid:
            return structure
    return None

def isNumber(s):
    try:
        float(s.replace(",",""))
        return True
    except ValueError:
        return False # Boolean

def isConcentraton(s):
    return isNumber(s) or isPercent(s) # Boolean

def isPercent(s):
    if (s[-3:] == "w/v" or s[-3:] == "v/v") and "%" in s:
        return True
    elif len(s) > 1 and s[-1] == "%":
        try:
            float(s[:-1])
            return True
        except ValueError:
            return False
    else:
        return False # Boolean

def isCompound(word): # boolean
    """Takes a string and returns true if the string is part of a chemical compound
    Returns false if the string is a comma or a number"""
    isPunctuation = word in [",", ".", ":", ";", "-"]
    return (not isConcentraton(word) and not isPunctuation and not word == "ERROR" and word not in STOP_WORDS)

def averageNumberString(numberRange): # float
    """Takes a string such as 20-25 and averages the two (or more) numbers (ex. "20-25" --> "22.5")
    Returns None if the split doesn't produce all numbers (ex. 20%-25%)
    """
    stringSplit = numberRange.split('-')
    if all(isNumber(number) for number in stringSplit):
        numbers = [float(i.replace(",","")) for i in stringSplit]
        return sum(numbers) / len(numbers)
    else:
        return None

def wordReplacement(s, replacementDictionary, debug=False):
    """Replaces certain strings to make parsing easier
    replacementDictionary: an OrderedDict dictionary which controls the replacements to make
    If debug is True, then every step is printed"""
    details = s
    ignoredSequence = "$*"
    for key, value in replacementDictionary.items():
        details = details.replace(key.replace(ignoredSequence, ""), value)
        if debug:
            print("Replacing '{}' with '{}':".format(key.replace(ignoredSequence, ""),value))
            print(details+"\n")
    return details

def loadPdb(pdbid): # ElementTree root
    """Loads a pdbid as an xml file and returns the <root> branch"""
    try:
        response = requests.get("http://www.rcsb.org/pdb/rest/customReport.xml?pdbids="+pdbid+
        "&customReportColumns=crystallizationMethod,crystallizationTempK,pdbxDetails,phValue,pmc,sequence,resolution&service=wsfile", timeout=10)
    except requests.Timeout:
        print("Request timeout on PDB: {}".format(pdbid))
        return None
    except KeyboardInterrupt:
        print("HTTP Request interrupted - exiting program")
        sys.exit()
    #sleep(0.1) # Limit HTTP requests
    root = etree.fromstring(response.content)
    if (root.find("record") != None):
        return root
    else:
        print("\n-------------------- ERROR --------------------\nNo entry found with PDB ID "
        + str(pdbid)+"\n-----------------------------------------------\n")
        return None

def fetchStructures(pdbList, structureFile, onlyDetails=True): # list
    """Takes a list of pdbids and creates Structure objects for them, outputing them to structureFile
    Also will read structureFile and ignore pdbs which have already been created
    If onlyDetails is true the function will only output Structures that have crystallization details
    """
    count = 1
    structureList = []

    pdbsWithoutDetails = []
    try:
        print("Loading PDBs without details...")
        pdbsWithoutDetails = loadJson(WITHOUT_DETAILS_FILE)
    except:
        print("File pdbs_without_details.json not found. One will be created.")

    completedPdbList = [] # List of Pdbs already turned into structures or without details
    try:
        structureList = loadStructures(structureFile)
        for struc in structureList:
            completedPdbList.append(struc.pdbid)
    except FileNotFoundError:
        print("File {} not found. A new file will be created.".format(structureFile))

    if onlyDetails:
        completedPdbList.extend(pdbsWithoutDetails)
    print("Removing completed pdbs...")
    pdbList = [pdbid for pdbid in pdbList if pdbid not in completedPdbList]

    for pdbid in pdbList:
        if count % 100 == 0 or count == 1:
            print("Loading pdb {} of {}...".format(count, len(pdbList)))
        root = loadPdb(pdbid)
        count += 1
        if (root != None):
            details = root.find("record/dimStructure.pdbxDetails").text
            if details == "null":
                details = None
                pdbsWithoutDetails.append(pdbid)
                writeJson(pdbsWithoutDetails, WITHOUT_DETAILS_FILE)
            if details != None or not onlyDetails:
                pmcid = root.find("record/dimStructure.pmc").text
                if (pmcid != "null"):
                    pmcid = pmcid[3:]
                else:
                    pmcid = None
                try:
                    pH = float(root.find("record/dimStructure.phValue").text)
                except ValueError:
                    pH = None
                try:
                    temperature = float(root.find("record/dimStructure.crystallizationTempK").text)
                except ValueError:
                    temperature = None
                method = root.find("record/dimStructure.crystallizationMethod").text
                if method == "null":
                    method = None
                sequences = []
                for tag in root.findall("record/dimEntity.sequence"):
                    sequences.append(tag.text)
                try:
                    resolution = float(root.find("record/dimStructure.resolution").text)
                except ValueError:
                    resolution = None
                structure = Structure(pdbid, pmcid, details, [], pH, temperature, method, sequences, resolution)
                structureList.append(structure)
                if count % 100 == 0:
                    writeStructures(structureList, structureFile)

    writeStructures(structureList, structureFile)
    print("Done fetching Structures")
    return structureList

def loadStructures(structureFile=STRUCTURES_FILE): # list
    """Returns a list of structures from the pickled structure file"""
    print("Loading structures from file...")
    return pickle.load(open(structureFile, "rb"))

def writeStructures(structureList, structureFile):
    """Writes a list of structures to a pickle file"""
    try:
        pickle.dump(structureList, open(structureFile, "wb"))
    except KeyboardInterrupt:
        print("Writing to file...")
        pickle.dump(structureList, open(structureFile, "wb"))
        sys.exit()

if __name__ == "__main__":
    structureList = loadStructures(STRUCTURES_FILE)
    # parseAllDetails(structureList, searchString="),")
    # exportSensibleStructures(structureList)
    # updateSmilesDictionary(structureList)
    # exportErrors(structureList, ERROR_DETAILS_FILE, unknown=True)
    # exportAllDetails(structureList, DETAILS_FILE)
    # exportCompoundFrequencies(structureList, COMPOUND_FREQUENCY_FILE)
    # exportCompoundFrequencies(structureList, UNKNOWN_FREQUENCY_FILE, mode="unknown")
    # exportCompoundFrequencies(structureList, PENDING_FREQUENCY_FILE, mode="pending")
    # getStructure(structureList, "1B2W").parseDetails(debug=True)
