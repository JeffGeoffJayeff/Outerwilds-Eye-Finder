# Program that lets the user type in commands to analyze data
# V 0.1

import re
import numpy as np
from pathlib import Path
from tabulate import tabulate #For making tables 
from plot_spherical import spherical_to_cartesian, cartesian_to_spherical
from probe import resultsDType

eyeShellRadius = 286500 #The radius of the eye shell in meters https://www.reddit.com/r/outerwilds/comments/t7mxcy/how_far_away_is_the_eye_base_game_spoilers/




class terminal:

    def __init__(self):
        self.dataset = [] #
        self.running = True
        self.eyeShellCartesian = None
        self.simulations = 0
        self.commands = {
            "quit": self.quit,
            "FileLoad": self.loadFile,
            "FolderLoad": self.loadFolder,
            "VisitStats": self.visitStats,
            "Help": self.helpCommand,
            "Print": self.printData,
            "Save": self.saveData,
            "Test": self.test,
            "AddCartesian": self.calculateCartesian,
            "SortbyVisits": self.sortbyVisits,
            "SumVisits": self.sumVisits,
            "LookupLaunch": self.lookupLaunchConditions,
            "LookupIndex": self.lookupIndex
        }
        self.resultsFields = resultsDType
        self.visitFields = [
        'Sun Visits',
        'Sun Station Visits',
        'Ember Twin Visits',
        'Ash Twin Visits',
        'Timber Hearth Visits',
        'Attlerock Visits',
        'Brittle Hollow Visits',
        "Hollow's Lantern Visits",
        "Giant's Deep Visits",
        'Cannon Visits',
        'Dark Bramble Visits',
        'Interloper Visits',
        'White Hole Visits',
        'Stranger Visits',
        'Random Eye Visits',
        'Spacey Visits'
        ] #Names of visit fields in teh dataset, I don't know if this changed in the various folders but they are correct for the latest two, UniformDistDifferentSpeed and UniformDistEyeHasMass
        self.totalVisitsName = "Total Visits" #name of the column that holds the sum of visits
    def commandRunner(self, input:str):
        splitInput = input.split(" ")
        command = splitInput[0]
        parameters = splitInput[1:]
        if command in self.commands:
            if len(parameters) == 0:
                self.commands[command]()
            elif len(parameters) > 0:
                self.commands[command](*parameters)
            else:
                print("ERROR: Invalid number of parameters")
        else:
            print("ERROR: Unknown command")
    def addColumn(self,newDataName:str,newDataType,newData = None): #Adds a column to
        self.resultsFields = self.resultsFields + [(newDataName,newDataType)]
        newArray = np.empty(len(self.dataset),dtype=self.resultsFields)
        for field in self.dataset.dtype.names:
            newArray[field] = self.dataset[field]
        if newData is None:
            self.dataset = newArray
            return #Datatype could be anything, and adding special initialization for the new column is not worth the effort, so just leave empty
        if np.size(newData) == len(self.dataset):
            newArray[newDataName] = newData
            self.dataset = newArray
        else:
            print(f"ERROR: New data length {np.size(newData)} does not match dataset length {len(self.dataset)}")
            return
        
    def calculateCartesian(self):
        shellx = np.cos(self.dataset['Eye Shell Polar']) * np.sin(self.dataset['Eye Shell Azimuth']) * eyeShellRadius
        shelly = np.sin(self.dataset['Eye Shell Polar']) * np.sin(self.dataset['Eye Shell Azimuth']) * eyeShellRadius
        shellz = np.cos(self.dataset['Eye Shell Azimuth']) * eyeShellRadius

        finalx = np.cos(self.dataset['Final Polar']) * np.sin(self.dataset['Final Azimuth']) * self.dataset['Final Radius']
        finaly = np.sin(self.dataset['Final Polar']) * np.sin(self.dataset['Final Azimuth']) * self.dataset['Final Radius']
        finalz = np.cos(self.dataset['Final Azimuth']) * self.dataset['Final Radius']

        self.addColumn("Eye Shell X",np.float64,shellx)
        self.addColumn("Eye Shell Y",np.float64,shelly)
        self.addColumn("Eye Shell Z",np.float64,shellz)

        self.addColumn("Final X",np.float64,finalx)
        self.addColumn("Final Y",np.float64,finaly)
        self.addColumn("Final Z",np.float64,finalz)
    def loadFile(self, filename):
        data = np.load(filename)
        self.dataset.append(data)

    def loadFolder(self, foldername): #TODO: Check that folder exists before doing this
        print(f"Loading folder {foldername}...",end="")
        npy_files = list(Path(foldername).glob("*.npy"))
        seperateData = []
        for filename in npy_files:
            seperateData.append(np.load(filename))
        combinedData = np.concatenate(seperateData)
        if len(self.dataset) > 0:
            self.dataset = np.concatenate([self.dataset, combinedData]) 
        else:
            self.dataset = combinedData
        self.simulations += np.size(combinedData,0)
        print(f"Folder {foldername} loaded with {len(seperateData):,d} files and {np.size(combinedData,0):,d} simulations\nTotal number of simulations: {self.simulations:,d}")
    def lookupIndex(self,index:int):
        index = int(index)
        if index < 0 or index >= len(self.dataset):
            print(f"ERROR: Index {index} is out of bounds for dataset of length {len(self.dataset)}")
            return
        row = self.dataset[index].tolist()
        table = [list(self.dataset.dtype.names),row]
        table = np.rot90(np.array(table))
        outputTable = tabulate(table, tablefmt="pretty")
        print(outputTable)
        #print(table)
    def lookupLaunchConditions(self,x,y,z,searchtype:str=None): #Look up the launch conditions of the simulation that results in the closest x y z coordinates on either the eye shell or the final probe position
        #Outputs as a numpy array of [unitx,unity,unitz,velocity,index]
        x = float(x)
        y = float(y)
        z = float(z)
        print("Checking/adding columns...",end="")
        if 'Distance to Point' not in self.dataset.dtype.names:
            self.addColumn('Distance to Point',np.float64) #This column stores the distance of every simulation to the input point
        if 'Eye Shell X' not in self.dataset.dtype.names:
            self.calculateCartesian() #Add the eye shell cartesian coordinates and final cartesian coordinates if they don't exist yet
        print("Checking search type...",end="")
        if searchtype.casefold() == "final".casefold():
            self.dataset['Distance to Point'] = np.sqrt((self.dataset['Final X'] - x)**2 + (self.dataset['Final Y'] - y)**2 + (self.dataset['Final Z'] - z)**2)
            self.dataset['Distance to Point'][np.isnan(self.dataset['Distance to Point'])] = np.inf #Set the distance of any calculations that resulted in NaN to infinity to avoid them showing up
            index = np.argmin(self.dataset['Distance to Point'])
        elif searchtype.casefold() == "eye".casefold():
            booleanmask = self.dataset['Reached Eye']
            #Set the distance to point of simulations that don't reach the eye to infinity so they are not considered
            self.dataset['Distance to Point'] = np.sqrt((self.dataset['Eye Shell X'] - x)**2 + (self.dataset['Eye Shell Y'] - y)**2 + (self.dataset['Eye Shell Z'] - z)**2)    
            self.dataset['Distance to Point'][np.logical_not(booleanmask)] = np.inf
            index = np.argmin(self.dataset['Distance to Point'])
        else: 
            print("ERROR: Invalid search type, must be either 'final' or 'eye'")
            return
        print("Outputting results...")
        unitx = self.dataset['Relative Launch x'][index]
        unity = self.dataset['Relative Launch y'][index]
        unitz = self.dataset['Relative Launch z'][index]
        velocity = self.dataset['Relative Launch Velocity'][index]
        print(f"Closest simulation to point ({x},{y},{z}) is at index {index} with the following parameters: \nDistance to point: {self.dataset['Distance to Point'][index]}\nLaunch velocity {velocity}\nUnit vector ({unitx}, {unity}, {unitz})")
        return np.array([unitx,unity,unitz,velocity,index])
    
    def helpCommand(self):
        print("Current Commands:")
        for cmd in self.commands:
            print(f"  {cmd}")

    def printData(self):
        print(self.dataset)
        print(type(self.dataset))


    def quit(self):
        self.running = False
        print("Bye bye!")
    
    def saveData(self, fileName:str, folderName:str=None):
        if len(self.dataset) > 0:
            if folderName is None:
                savePath = f"{fileName}.npy"
            else:
                savePath = f"{folderName}/{fileName}.npy"
            np.save(savePath, self.dataset)
            print(f"Data saved to {savePath}")
        else:
            print("No data to save")
    def sortbyVisits(self):
        if any(self.totalVisitsName != name for name in self.dataset.dtype.names):
            print(f"{self.totalVisitsName} not found, calculating sum of visits...")
            self.sumVisits()
        print("Sorting dataset by total visits in descending order...")
        self.dataset = self.dataset[np.argsort(self.dataset[self.totalVisitsName])[::-1]]
        print("Dataset sorted by total visits in descending order")
    def sumVisits(self):
        visit_fields = [name for name in self.dataset.dtype.names if 'Visits' in name]
        sumColumn = np.sum([self.dataset[field] for field in visit_fields], axis=0)
        self.addColumn(self.totalVisitsName, np.int32, sumColumn)
    def test(self):
        return

    def visitStats(self):
        visitnums = []
        totalVisits = 0

        # Calculating hit values
        unique, counts = np.unique(self.dataset['Body Hit'], return_counts=True)
        totalHits = np.sum(self.dataset['Hit Something'])

        countDict = dict(zip(unique, counts))
        for body in self.visitFields:
            totalVisits += np.sum(self.dataset[body]) #Bit wasteful to sum twice but it is clearer
        for index, body in enumerate(self.visitFields):
            
            visitSum = np.sum(self.dataset[body])
            visitPercent = visitSum / totalVisits * 100 #Percentage of visits to this body out of all visits to all bodies
            visitofAllSims = visitSum / self.simulations * 100 #Percentage of visits to this body out of all simulations

            if index in countDict:
                hitSum = countDict[index] #Doing this because np.unique doesn't return the index of bodies that were never hit
            else:
                hitSum = 0
            hitPercent = hitSum / totalHits * 100 #Percentage of hits to this body out of all hits to all bodies
            hitofAllSims = hitSum / self.simulations * 100 #Hits to this body of all simulations

            tableEntry = [re.sub(" Visits","",body)] #Remove visits from the names of the bodies for the table
            tableEntry.append(f"{visitSum:,}")
            if visitSum != 0:
                tableEntry.append(f"{visitPercent:.3f}")
                tableEntry.append(f"{visitofAllSims:.3f}")
            else:
                tableEntry.append("0")
                tableEntry.append("0")

            tableEntry.append(f"{hitSum:,}") #Add the number of hits to the table
            if hitSum != 0:
                tableEntry.append(f"{hitPercent:.3f}")
                tableEntry.append(f"{hitofAllSims:.3f}")
            else:
                tableEntry.append("0")
                tableEntry.append("0")

            visitnums.append(tableEntry)
        outputTable = tabulate(visitnums,showindex=True,headers=["Body","Visits","Visit/Visits %","Visit/Sim %","Hits","Hits/Hits %","Hit/Sim %"],tablefmt="pretty")
        print(outputTable)
def main():
    print("Start of session")
    termGuy = terminal()
    while termGuy.running:
        termGuy.commandRunner("Help")
        command = input("Enter a command (or 'quit' to exit): ")
        termGuy.commandRunner(command)

main()