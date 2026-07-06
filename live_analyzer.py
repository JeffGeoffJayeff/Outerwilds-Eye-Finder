# Program that lets the user type in commands to analyze data
# V 0.1

import re
import numpy as np
from pathlib import Path
from tabulate import tabulate #For making tables 
from plot_spherical import spherical_to_cartesian, cartesian_to_spherical


class terminal:

    def __init__(self):
        self.dataset = [] #
        self.running = True
        self.simulations = 0
        self.commands = {
            "quit": self.quit,
            "FileLoad": self.loadFile,
            "FolderLoad": self.loadFolder,
            "VisitStats": self.visitStats,
            "Help": self.helpCommand,
            "Print": self.printData,
            "Save": self.saveData
        }
        self.visit_fields = [
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

    def loadFile(self, filename):
        data = np.load(filename)
        self.dataset.append(data)

    def loadFolder(self, foldername):
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

    def lookupLaunchConditions(self, x,y,z): #Look up the launch conditions of the simulation that results in the closest x y z coordinates on the eye shell
        desiredCoords = cartesian_to_spherical(x,y,z)
        
    def visitStats(self):
        visitnums = []
        totalVisits = 0

        # Calculating hit values
        unique, counts = np.unique(self.dataset['Body Hit'], return_counts=True)
        totalHits = np.sum(self.dataset['Hit Something'])

        countDict = dict(zip(unique, counts))
        for body in self.visit_fields:
            totalVisits += np.sum(self.dataset[body]) #Bit wasteful to sum twice but it is clearer
        for index, body in enumerate(self.visit_fields):
            
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
    def helpCommand(self):
        print("Current Commands:")
        for cmd in self.commands:
            print(f"  {cmd}")

    def printData(self):
        print(self.dataset)
        
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
def main():
    print("Start of session")
    termGuy = terminal()
    while termGuy.running:
        termGuy.commandRunner("Help")
        command = input("Enter a command (or 'quit' to exit): ")
        termGuy.commandRunner(command)

main()