import pickle, json, os, sys
import xml.etree.ElementTree as etree
from pathlib import Path
from pdb_crystal_database import loadStructures, writeStructures, getStructure, Structure
from misc_functions import loadJson, writeJson

try:
    import requests
except ModuleNotFoundError:
    print("ERROR: Requests module not found. You must install the requests module in order to get structure info from the PDB."
    "For more information, see http://docs.python-requests.org/en/master/user/install/")
    sys.exit()

if not os.path.exists("Structures"):
    os.makedirs("Structures")

# Create Path objects for directories
INPUT_DIR = Path("Input/")
STRUCTURE_DIR = Path("Structures/")

WITHOUT_DETAILS_FILE = INPUT_DIR / "pdbs_without_details.json"
STRUCTURES_FILE = STRUCTURE_DIR / "structures.pkl" # The database file. Must be placed in proper location

def loadPdb(pdbid): # ElementTree root
    """Loads a pdbid as an xml file and returns the <root> branch"""
    try:
        response = requests.get("http://www.rcsb.org/pdb/rest/customReport.xml?pdbids="+pdbid+
        "&customReportColumns=crystallizationMethod,crystallizationTempK,pdbxDetails,phValue,pmc,sequence,resolution&service=wsfile", timeout=10)
    except (requests.Timeout, requests.exceptions.ConnectionError):
        print("Request timeout on PDB: {}".format(pdbid))
        return None
    except KeyboardInterrupt:
        print("HTTP Request interrupted - exiting program")
        sys.exit()
    root = etree.fromstring(response.content)
    if (root.find("record") != None):
        return root
    else:
        print("\n-------------------- ERROR --------------------\nNo entry found with PDB ID "
        + str(pdbid)+"\n-----------------------------------------------\n")
        return None

def fetchStructures(pdbList, structureFile=STRUCTURES_FILE, onlyDetails=True, ignorePdbsWithoutDetails=True, ignoreCompletedPdbs=True, saveFrequency=100): # list
    """Takes a list of pdbids and creates Structure objects for them, outputing them to structureFile
    If onlyDetails is True the function will only output Structures that have crystallization details
    If ignoreCompletedPdbs is True, then the function will read the Structure file and ignore Pdbs which already
        have Structures associated with them
    If ignoreCompletedPdbs is False, then the function will recheck and update every structure in the structure file
    If ignorePdbsWithoutDetails is True then the function will ignore Pdbs from WITHOUT_DETAILS_FILE
    saveFrequency is the number of pdbs downloaded before the files are saved
        increasing this number makes the script faster, but more unsaved data is lost when the script is exited
    """
    count = 1
    structureList = []

    pdbsWithoutDetails = []
    try:
        print("Loading PDBs without details...")
        pdbsWithoutDetails = loadJson(WITHOUT_DETAILS_FILE)
    except:
        print("File {} not found. One will be created.".format(WITHOUT_DETAILS_FILE))

    completedPdbList = [] # List of Pdbs already turned into structures or without details
    try:
        structureList = loadStructures(structureFile)
        if ignoreCompletedPdbs:
            for struc in structureList:
                completedPdbList.append(struc.pdbid)
    except FileNotFoundError:
        print("File {} not found. A new structure file will be created.".format(structureFile))

    ignoredPdbList = []

    if ignorePdbsWithoutDetails:
        ignoredPdbList.extend(pdbsWithoutDetails)
    # If ignoreCompletedPdbs=False, then completedPdbList will be empty
    ignoredPdbList.extend(completedPdbList)

    print("Removing completed pdbs (May take a minute)...")
    pdbList = list(set(pdbList)-set(ignoredPdbList))

    if pdbList == []:
        print("All PDBs have been ignored. No more PDBs need to be downloaded.")

    for pdbid in pdbList:
        if count == 1:
            print("Downloading {} structure objects from the pdb".format(len(pdbList)))
        if count % 100 == 0:
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
                # If the pdb already has a structure in the list, update it
                structureAlreadyInList = getStructure(structureList, structure.pdbid)
                if structureAlreadyInList != None:
                    structureList.remove(structureAlreadyInList)
                structureList.append(structure)
                if count % saveFrequency == 0:
                    # Save structures
                    writeStructures(structureList, structureFile)

    writeStructures(structureList, structureFile)
    print("Done fetching Structures")
    return structureList

def getAllPdbs(filename=""): # list
    """Returns a list of every pdbid in the PDB
    # filename = an optional argument to also export the list as a json file"""
    print("Getting a list of all pdbs...")
    response = requests.get("https://www.rcsb.org/pdb/json/getCurrent").text
    dictionary = json.loads(response)
    allPdbs = dictionary["idList"]
    if filename != "":
        writeJson(allPdbs, filename)
    return allPdbs

if __name__ == "__main__":
    allPdbs = getAllPdbs()
    fetchStructures(allPdbs)
