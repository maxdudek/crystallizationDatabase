import json, csv, sys, pickle, collections
from nltk import word_tokenize
from misc_functions import loadJson, writeJson, printList, getKey, listToFile, fileToList
from pdb_crystal_database import Structure
from pdb_crystal_database import loadStructures, parseAllDetails, writeStructures

UNKNOWN_LIST_FILE = "Input\\unknown_list.json"
COMPOUND_DICTIONARY_FILE = "Input\\compound_dictionary.json"
STRUCTURES_FILE = "Structures\\structures.pkl"
STOP_WORDS_FILE = "Input\\stop_words.json"
ERROR_LIST_FILE = "Input\\error_list.json"

print("Loading dictionary files...")
compoundDictionary = loadJson(COMPOUND_DICTIONARY_FILE)
stopWords = loadJson(STOP_WORDS_FILE)
unknownList = loadJson(UNKNOWN_LIST_FILE)
errorList = loadJson(ERROR_LIST_FILE)
structureList = loadStructures(STRUCTURES_FILE)

# Define input options
INPUT_QUIT = "" # Exit the script
INPUT_SAVE = "save" # Save files - currently done automatically after every change
INPUT_UNKNOWN = "unknown" # Add the compound to the unknownList
INPUT_SAME = "=" # Add the compound to the dictionary exactly as it appears (e.g. "sodium chloride" --> "sodium chloride")
INPUT_UNDO = "u" # Undo - NOTE: undo is unreliable and should only be used for single dictionary changes
# To undo a mistake safely, EXIT THE SCRIPT and make the change in the text files directly
INPUT_IGNORE = "sw" # Add ALL WORDS in the compound to the list of stop words
# ex. if the compound is "well plate", both "well" and "plate" are added to the list
INPUT_ADD_STOP_WORD = "add" # Adds a single stop word to the stopWord list
# ex. "add well" will append "well" to the stopWords list
INPUT_PASS = "pass" # Skips the current compound, does nothing to it
INPUT_REMOVE_STOP_WORDS = "rm" # removes stop words from the compound and reloads the prompt
INPUT_ERROR = "error" # Adds the compound to the error list


def getCompoundList(structureList, sortedByFrequency=True, getKey=False): # list
    """Takes in a list of Structure Objects and returns a list of just the compound names
    If sortedByFrequency is True, the compound names are sorted by frequency
    If getKey is true, then the compounds are turned into dictionary keys using getKey(), and then sorted (eg "na-cl" --> "nacl")
    """
    print("Getting compound list...")
    compoundList = []
    for structure in structureList:
        compoundList.extend(structure.compounds[0::2])
    if getKey:
        compoundList = [getKey(compound) for compound in compoundList if compoundList]
    counts = collections.Counter(compoundList)
    if sortedByFrequency:
        compoundList = sorted(compoundList, key=lambda x: -counts[x])
    return compoundList

def removeStopWords(s, stopWords):
    """Takes a string and returns a string with stop words removed
    stopWords = A list of words to remove"""
    words = word_tokenize(s)
    words = [word for word in words if word not in stopWords]
    return printList(words, " ")

def saveFiles():
    """Saves all of the lists and dictionaries to their respective files"""
    writeJson(compoundDictionary, COMPOUND_DICTIONARY_FILE, indent=2)
    writeJson(unknownList, UNKNOWN_LIST_FILE, indent=2)
    writeJson(errorList, ERROR_LIST_FILE, indent=2)
    writeJson(stopWords, STOP_WORDS_FILE, indent=2)
    print("Files saved")

def generateDictionary(compoundList): # dictionary
    """ Iterates through a list and substitutes elements based on a dictionary
    If no key is found for the element, the user is prompted to enter an entry
    See INPUT definitions above for more options
    """
    print("Beginning dictionary generation (may take a minute)...")
    history = [] # A list of modified indeces, in order to UNDO
    i = 0
    while(i < len(compoundList)+1):
        if i < len(compoundList):
            compound = compoundList[i]
            compound = removeStopWords(compound, stopWords)
            if getKey(compound) in compoundDictionary or getKey(compound) in unknownList or getKey(compound) in errorList:
                pass
            elif compound in [" ", "", "-", ":"]:
                pass
            else: # Name of compound not found in dictionary or lists
                print("Reading compound {} of {}".format(i+1, len(compoundList)))
                runAgain = True
                while(runAgain and getKey(compound) not in compoundDictionary):
                    runAgain = False
                    inputText = input("Enter the name of the following compound:\n{}\n$:".format(compound))
                    # PARSE INPUT
                    if inputText == INPUT_QUIT: # Quit
                        inputText = input("Really quit?\n0: Save and quit\n1: Quit without saving\n2: Don't quit\n$:")
                        if inputText == "0":
                            saveFiles()
                            sys.exit()
                        elif inputText == "1":
                            sys.exit()
                        else:
                            runAgain = True
                    elif inputText == INPUT_SAVE: # Save
                        saveFiles()
                        runAgain = True
                    elif inputText == INPUT_UNKNOWN: # Unknown compound
                        unknownList.append(getKey(compound))
                        history.append(i)
                    elif inputText == INPUT_IGNORE:
                        ignored_words = word_tokenize(compound)
                        for word in ignored_words:
                            if word not in stopWords:
                                stopWords.append(word)
                                writeJson(stopWords, STOP_WORDS_FILE, indent=2)
                        history.append(i)
                    elif inputText == INPUT_ERROR:
                        errorList.append(getKey(compound))
                        history.append(i)
                    elif inputText == INPUT_UNDO:
                        if history == []:
                            print("Unable to undo")
                            runAgain = True
                        else:
                            oldIndex = history[-1]
                            oldNameKey = getKey(compoundList[oldIndex])
                            if oldNameKey in compoundDictionary:
                                oldValue = compoundDictionary[oldNameKey]
                                del compoundDictionary[oldNameKey]
                                print("Removed key from dictionary:\n{} : {}".format(oldNameKey, oldValue))
                            elif oldNameKey in unknownList:
                                del unknownList[index(oldNameKey)]
                                print("Removed {} from unknownList".format(oldNameKey))
                            elif oldNameKey in errorList:
                                del errorList[index(oldNameKey)]
                                print("Removed {} from errorList".format(oldNameKey))
                            saveFiles()
                            del history[-1]
                            i = oldIndex - 1
                    elif inputText[:len(INPUT_ADD_STOP_WORD)] == INPUT_ADD_STOP_WORD: # Add stop word
                        wordToAdd = inputText[len(INPUT_ADD_STOP_WORD)+1:]
                        if wordToAdd not in stopWords:
                            stopWords.append(wordToAdd)
                            writeJson(stopWords, STOP_WORDS_FILE, indent=2)
                            print("Added stop word {}".format(wordToAdd))
                            compound = removeStopWords(compound, stopWords)
                            runAgain = True
                    elif inputText == INPUT_REMOVE_STOP_WORDS:
                        compound = removeStopWords(compound, stopWords)
                        runAgain = True
                    elif inputText == INPUT_PASS:
                        pass
                    else: # Normal input
                        nameOfCompound = inputText
                        if inputText == INPUT_SAME:
                            nameOfCompound = compound
                        inputText = input("Add the following key to the dictionary? (Press n to cancel, ENTER to confirm):\n{} : {}\n$:".format(compound, nameOfCompound))
                        if inputText != "n":
                            compoundDictionary[getKey(compound)] = nameOfCompound
                            # Add the value to the dictionary with itself as the key
                            if getKey(nameOfCompound.lower()) not in compoundDictionary:
                                compoundDictionary[getKey(nameOfCompound.lower())] = nameOfCompound
                            history.append(i)
                            print("Added")
                        else:
                            runAgain = True
                saveFiles()
        else: # END
            runAgain = True
            while(runAgain):
                runAgain = False
                inputText = input("Reached end of list. Save to dictionary? (Press n to cancel, ENTER to confirm, u to undo)\n$:")
                if inputText == "n":
                    inputText = input("Are you sure you want to quit without saving? (y/n)\n$:")
                    if inputText == "y":
                        sys.exit()
                    else:
                        saveFiles()
                elif inputText == INPUT_UNDO: # Undo
                    if history == []:
                        print("Unable to undo")
                        runAgain = True
                    else:
                        oldIndex = history[-1]
                        oldNameKey = getKey(compoundList[oldIndex])
                        if oldNameKey in compoundDictionary:
                            oldValue = compoundDictionary[oldNameKey]
                            del compoundDictionary[oldNameKey]
                            print("Removed key from dictionary:\n{} : {}".format(oldNameKey, oldValue))
                        elif oldNameKey in unknownList:
                            del unknownList[index(oldNameKey)]
                            print("Removed {} from unknownList".format(oldNameKey))
                        elif oldNameKey in errorList:
                            del errorList[index(oldNameKey)]
                            print("Removed {} from errorList".format(oldNameKey))
                        saveFiles()
                        del history[-1]
                        i = oldIndex - 1
                else: # Any other input
                    saveFiles()
                    sys.exit()
        i += 1


# parseAllDetails(structureList)
# writeStructures(structureList, STRUCTURES_FILE)
compoundList = getCompoundList(structureList)
try:
    generateDictionary(compoundList)
except KeyboardInterrupt:
    saveFiles()
    sys.exit()
