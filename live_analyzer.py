# Program that lets the user type in commands to analyze data
# V 0.1

import numpy as np
from pathlib import Path

class terminal:

    def __init__(self):
        self.dataset = [] #
        self.running = True
        self.commands = {
            "quit": self.quit,
            "FileLoad": self.loadFile,
            "FolderLoad": self.loadFolder,
            "VisitStats": self.visitStats,
            "help": self.helpCommand,
            "print": self.printData
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
            else:
                self.commands[command](*parameters)
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
        self.dataset = combinedData
        print(f"Folder {foldername} loaded with {len(seperateData):,d} files and {np.size(self.dataset,0):,d} simulations")

    def visitStats(self):
        vistnums = []
        for body in self.visit_fields:
            vistnums.append(np.sum(self.dataset[body]))
        print(vistnums)

    def helpCommand(self):
        print("Current Commands:")
        for cmd in self.commands:
            print(f"  {cmd}")

    def printData(self):
        print(self.dataset)
        
    def quit(self):
        self.running = False
        print("Bye bye!")
def main():
    print("Start of session")
    termGuy = terminal()
    while termGuy.running:
        termGuy.commandRunner("help")
        command = input("Enter a command (or 'quit' to exit): ")
        termGuy.commandRunner(command)

main()