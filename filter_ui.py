import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showerror, showinfo, askyesno
import datetime
import webbrowser
import gzip
import pathlib
from typing import Union, Any, List
from multiprocessing.pool import ThreadPool as Pool
import threading

import matplotlib as plt
plt.use('TkAgg')

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)

Path = Union[str, pathlib.Path]

class ThreadedParser():
    def __init__(self, path: Path): 
        self.path = path

    def getNextLine(self):
        f = open(self.path)
        
        if self.path.endswith('.gz'):
            f = gzip.open(self.path, 'rt') 

        for line in f:
            yield line

    def processLine(self, line: str): 
        pass

    def afterProcess(self) -> Any: 
        pass

    def run(self):
        pool = Pool(processes=8)

        for line in self.getNextLine():
            pool.map(self.processLine, (line, ))
        
        pool.close()
        pool.join()

        return self.afterProcess()

class MetaDataParser(ThreadedParser):
    def __init__(self, path: Path, phenotypes: List[str]):
        super().__init__(path)
        self.found = dict.fromkeys(phenotypes)
        for k in self.found.keys():
            self.found[k] = []

    def processLine(self, line: str):
        data = tuple(x.strip() for x in line.split("\t"))

        for phenotype in self.found.keys():
            if phenotype.lower() in line.lower():
                self.found[phenotype].append((data[0], data[3]))
                
    def afterProcess(self):
        return self.found

class ChrParser(ThreadedParser):
    def __init__(self, path: Path, chr: str, pos: str, ref: str, alt: str):
        super().__init__(path)

        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt

        self.found = {}

    def processLine(self, line):
        data = tuple(x.strip() for x in line.split('\t'))

        if data[0] != self.chr:
            return

        if data[1] != self.pos:
            return

        if data[2] != self.ref:
            return

        found = False    
        for alt in data[3].split(','):
            if alt == self.alt:
                found = True
                break

        if not found:
            return

        alts = tuple(x for x in data[3].split(',') if x)
        zygosities = tuple(x for x in data[4].split(';') if x)
        ids = tuple(x for x in data[5].split(';') if x)
        
        # 0/1: Het for first alt
        # 1/1: Hom for first alt

        # 0/2: Het for second alt
        # 2/2: Hom for second alt

        # 0/3: Het for third alt
        # 3/3: Hom for third alt

        for zygosity, id in zip(zygosities, ids): 
            self.found[id] = 'het' if '0' in zygosity else 'hom'

    def afterProcess(self): 
        return self.found

class AsyncSearch(threading.Thread):
    def __init__(self, chr, pos, ref, alt, phenotypes, chrFile, metaFile, progressCallback, finishCallback):
        super().__init__()
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.phenotypes = phenotypes
        self.chrFile = chrFile
        self.metaFile = metaFile
        self.progressCallback = progressCallback
        self.finishCallback = finishCallback

    
    def run(self):
        phenotypes = MetaDataParser(self.metaFile, self.phenotypes.split(',')).run()
        self.progressCallback(33)
        
        chromosomes = ChrParser(self.chrFile, self.chr, self.pos, self.ref, self.alt).run()
        self.progressCallback(66)

        result = {}

        result['phenotypes'] = phenotypes
        result['chromosomes'] = chromosomes
        

        for index, (phenotype, pids) in enumerate(phenotypes.items()): 
            hets = []
            homs = []
            other_phenotypes = {}

            for pid in pids:
                if pid[0] in chromosomes.keys():
                    if chromosomes[pid[0]] == 'het':
                        hets.append(pid[0])
                    else:
                        homs.append(pid[0])

                    other_phenotypes[pid[0]] = pid[1]
                    
            
            if len(pids) > 0:
                freq = (1 * len(hets) + 2 * len(homs)) / (2 * len(pids))
            else: 
                freq = 0

            
            result[phenotype] = (freq, len(hets), len(homs), hets, homs, other_phenotypes)

        self.progressCallback(100)
        self.finishCallback(result)

class LabeledEntry:
    def __init__(self, master, label, row): 
        self.label = ttk.Label(master, text=label)
        self.label.grid(column=0, row=row, pady=2)
        self.entry = ttk.Entry(master)
        self.entry.grid(column=1, row=row, sticky='we', pady=2, padx=5, ipady=2, ipadx=2, columnspan=2)

class ReadonlyEntry:
    def __init__(self, master, label, row, textvariable = None): 
        self.label = ttk.Label(master, text=label)
        self.label.grid(column=0, row=row, pady=2)
        self.entry = ttk.Entry(master, state="readonly", textvariable=textvariable)
        self.entry.grid(column=1, row=row, sticky='we', pady=2, padx=5, ipady=2, ipadx=2, columnspan=2)

class MainWindow(tk.Tk): 
    def openChr(self): 
        filename = fd.askopenfilename(title='Open the meta data file', filetypes=(('Uncompressed txt', '*.txt'), ('Compressed txt', '*.gz')))

        if filename: 
            self.chrPath.set(filename)

            if self.metaDataPath.get() != 'No file is opened': 
                self.enableAfterOpen()

    def openMetaData(self):
        filename = fd.askopenfilename(title='Open the meta data file', filetypes=(('Uncompressed txt', '*.txt'), ('Compressed txt', '*.gz')))

        if filename: 
            self.metaDataPath.set(filename)

            if self.chrPath.get() != 'No file is opened': 
                self.enableAfterOpen()

    def enableAfterOpen(self):
        self.searchButton['state'] = 'normal'
        self.progress.stop()
        self.progress['mode'] = 'determinate'
        self.progress['value'] = 0

    def about(self):
        showinfo(title='About', message=f'Filteration UI\nMade by Reyhaneh Ahani\nCredit 2022-{datetime.date.today().year}')
    
    def progressCallback(self, value):
        pass

    def finishCallback(self, value):
        self.result = value
        self.phenotypes_result = value.pop('phenotypes')
        self.chromosomes_result = value.pop('chromosomes')

        self.searchButton['state'] = 'normal'
        self.openChrButton['state'] = 'normal'
        self.openMetaDataButton['state'] = 'normal'
        self.fileMenu.entryconfigure("Open Meta data", state='normal')
        self.fileMenu.entryconfigure("Open Chromosome file", state='normal')
        self.fileMenu.entryconfigure("Save plot as png", state='normal')
        self.progress['mode'] = 'determinate'
        self.progress.stop()
    
        self.chr.entry['state'] = 'normal'
        self.pos.entry['state'] = 'normal'
        self.ref.entry['state'] = 'normal'
        self.alt.entry['state'] = 'normal'
        self.phenotypes.entry['state'] = 'normal'

        self.savePngButton['state'] = 'normal'

        self.figure.clf()
        self.axes = self.figure.add_subplot()

        self.annot = self.axes.annotate("", xy=(0,0), xytext=(0,0),textcoords="offset points",
                                        bbox=dict(boxstyle="round", fc="white", ec="b", lw=2),
                                        arrowprops=dict(arrowstyle="->"))
        
        max_freq = 0
        for index, (key, items) in enumerate(self.result.items()): 
            freq = items[0]
            

            if freq > max_freq:
                max_freq = freq

            self.axes.bar(index, freq, width=1, edgecolor="white", linewidth=0.7)

        self.axes.set_ylim((0, max_freq + 0.1))
        self.axes.set_xlabel('Phenotypes')
        self.axes.set_ylabel('Frequency')
        self.axes.set_xticks(range(len(self.result.keys())), labels=self.result.keys(), rotation=45)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def savePlotAsPng(self):
        filename = fd.asksaveasfilename(title="Save plot", filetypes=(('PNG files', '*.png')))

        if filename:
            self.figure.savefig(filename)

    def monitor(self, thread): 
        if thread.is_alive():
            self.after(1000, lambda: self.monitor(thread))
        else:
            showinfo(title='Info', message=f'Filteration completed!')


    def processFiles(self): 
        if not (self.chr.entry.get() and 
                self.pos.entry.get() and
                self.ref.entry.get() and
                self.alt.entry.get() and
                self.phenotypes.entry.get()):
            showerror(title='Incomplete info', message='Please fill all required fields.')
        

        answer = askyesno(title='Confirmation', message='Searching in files is an expensive operation.\nAre you sure to proceed?')

        if not answer:
            return
        
        self.searchButton['state'] = 'disabled'
        self.openChrButton['state'] = 'disabled'
        self.openMetaDataButton['state'] = 'disabled'
        self.savePngButton['state'] = 'disabled'
        self.fileMenu.entryconfigure("Open Meta data", state='disabled')
        self.fileMenu.entryconfigure("Open Chromosome file", state='disabled')
        self.fileMenu.entryconfigure("Save plot as png", state='disabled')
        self.progress['mode'] = 'indeterminate'
        self.progress.start()
        self.chr.entry['state'] = 'disabled'
        self.pos.entry['state'] = 'disabled'
        self.ref.entry['state'] = 'disabled'
        self.alt.entry['state'] = 'disabled'
        self.phenotypes.entry['state'] = 'disabled'
        
        
        thread = AsyncSearch(self.chr.entry.get(),
                            self.pos.entry.get(),
                            self.ref.entry.get(),
                            self.alt.entry.get(),
                            self.phenotypes.entry.get(),
                            self.chrPath.get(),
                            self.metaDataPath.get(),
                            self.progressCallback,
                            self.finishCallback)

        thread.setDaemon(True)
        thread.start()
        self.monitor(thread)

    def plotLeave(self, event):
        pass
        
    def plotEnter(self, event):
        pass

    def plotHover(self, event):
        axes = event.inaxes
        
        if not hasattr(self, 'result'):
            return
        
        if axes:
            x = event.xdata  
            labels = axes.get_xticklabels()

            for index, patch in enumerate(axes.patches):
                if patch.get_bbox().containsx(x):
                    items = self.result[labels[index]._text]

                    freq = items[0]
                    nHet = items[1]
                    nHom = items[2]
                    hets = items[3]
                    homs = items[4] 

                    text = f'''Freq: {round(freq, 3)}\nn(Het): {nHet}\nn(Hom): {nHom}'''
                    if nHet > 0 or nHom > 0:
                        text += f'''\nHet IDS: {','.join(hets)}\nHom IDS: {','.join(homs)}'''
                        
                    self.annot.set_visible(True)
                    self.annot.set_text(text)
                    self.annot.xy = (index - 0.2, 0.05 + freq)
                    self.figure.canvas.draw_idle()

    def plotClick(self, event):
        axes = event.inaxes
        
        if not hasattr(self, 'result'):
            return
        
        if axes:
            x = event.xdata  
            labels = axes.get_xticklabels()

            for item in self.infoTable.get_children():
                self.infoTable.delete(item)

            for index, patch in enumerate(axes.patches):
                if patch.get_bbox().containsx(x):
                    items = self.result[labels[index]._text]
                    hets = items[3] 
                    for index, (key, item) in enumerate(items[5].items()):
                        if key in hets:
                            zygosity = 'Het'
                        else:
                            zygosity = 'Hom'
                        self.infoTable.insert('', 'end',text=str(index),values=(key, item, zygosity))

                    
    def __init__(self):
        super().__init__()

        self.metaDataPath = tk.StringVar(value='No file is opened')
        self.chrPath = tk.StringVar(value='No file is opened')

        self.title('Human Gene Pars :: Meta data filter UI')
        self.geometry('800x350')

        self.menuBar = tk.Menu(self)
        self.config(menu=self.menuBar)

        self.fileMenu = tk.Menu(self.menuBar, tearoff=False)
        self.fileMenu.add_command(label="Open Meta data", command=self.openMetaData)
        self.fileMenu.add_command(label="Open Chromosome file", command=self.openChr)
        self.fileMenu.add_command(label="Save plot as png", command=self.savePlotAsPng, state="disabled")
        self.fileMenu.add_separator()
        self.fileMenu.add_command(label="Exit", command=self.destroy)
        self.menuBar.add_cascade(label="File", menu=self.fileMenu)

        self.aboutMenu = tk.Menu(self.menuBar, tearoff=False)
        self.aboutMenu.add_command(label="Help", command=lambda: webbrowser.open('https://humangene.ir/'))
        self.aboutMenu.add_command(label="About", command=self.about)
        self.menuBar.add_cascade(label="About", menu=self.aboutMenu)

        self.fileManagement = ttk.LabelFrame(self, text='Options')
        self.fileManagement.pack(ipadx=10, ipady=10, fill='both', side='left', expand=True)

        self.fileManagement.columnconfigure(tuple(range(3)), weight=1)

        self.chr = LabeledEntry(self.fileManagement, 'Chromosome', 0)
        self.pos = LabeledEntry(self.fileManagement, 'Position', 1)
        self.ref = LabeledEntry(self.fileManagement, 'Reference', 2)
        self.alt = LabeledEntry(self.fileManagement, 'Alternative', 3)
        self.phenotypes = LabeledEntry(self.fileManagement, 'Phenotypes (seperated by comma)', 4)

        ttk.Separator(self.fileManagement, orient='horizontal').grid(column=0, row=5, columnspan=3, pady=5)
        
        self.openMetaDataButton = ttk.Button(self.fileManagement, text='Open meta data file', command=self.openMetaData)
        self.openMetaDataButton.grid(column=0, row=6, sticky='we', ipady=3, ipadx=3)

        self.openChrButton = ttk.Button(self.fileManagement, text='Open chromosome file', command=self.openChr)
        self.openChrButton.grid(column=1, row=6, sticky='we', ipady=3, ipadx=3)

        self.savePngButton = ttk.Button(self.fileManagement, text='Save plot as png', command=self.savePlotAsPng, state='disabled')
        self.savePngButton.grid(column=2, row=6, sticky='we', ipady=3, ipadx=3)

        ttk.Separator(self.fileManagement, orient='horizontal').grid(column=0, row=7, columnspan=3, pady=5)

        self.metaDataPathLabel = ReadonlyEntry(self.fileManagement, 'Meta data file:', 8, self.metaDataPath)
        self.chrPathLabel = ReadonlyEntry(self.fileManagement, 'Chromosome file:', 9, self.chrPath)
        
        ttk.Separator(self.fileManagement, orient='horizontal').grid(column=0, row=10, columnspan=3, pady=5)

        self.searchButton = ttk.Button(self.fileManagement, text='Search', state='disabled', command=self.processFiles)
        self.searchButton.grid(column=0, row=11, sticky='we', columnspan=1, ipady=3, ipadx=3)

        self.progress = ttk.Progressbar(self.fileManagement, orient='horizontal', length=100)
        self.progress.grid(column=1, row=11, sticky='we', columnspan=2, ipady=3, ipadx=3)

        self.infoTable = ttk.Treeview(self.fileManagement, column=('ID', 'Phenotypes', 'Zygosity'), show='headings', height=6)
        self.infoTable.column("# 1", anchor='center', width=75, stretch=False)
        self.infoTable.heading("# 1", text= 'ID')
        self.infoTable.column("# 2", anchor='center')
        self.infoTable.heading("# 2", text='Phenotypes')
        self.infoTable.column("# 3", anchor='center', width=75, stretch=False)
        self.infoTable.heading("# 3", text='Zygosity')  
        self.infoTable.grid(column=0, row=12, sticky='we', columnspan=3, ipady=3, ipadx=3)



        self.plotFrame = ttk.LabelFrame(self, text='Plot')
        self.plotFrame.pack(ipadx=10, ipady=10, fill='both', side='right', expand=True)
        self.figure = Figure(figsize=(6, 4), dpi=100)
        
        self.figure.canvas.mpl_connect('motion_notify_event', self.plotHover)
        self.figure.canvas.mpl_connect('axes_leave_event', self.plotLeave)
        self.figure.canvas.mpl_connect('axes_enter_event', self.plotEnter)
        self.figure.canvas.mpl_connect('button_press_event', self.plotClick)

        self.figureCanvas = FigureCanvasTkAgg(self.figure, self.plotFrame)
        NavigationToolbar2Tk(self.figureCanvas, self.plotFrame)
        self.figureCanvas.get_tk_widget().pack(side='top', fill='both', expand=True)


if __name__ == "__main__":
    app = MainWindow()
    app.mainloop()
