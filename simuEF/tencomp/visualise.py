import fedoo as fd

result = fd.read_data("simuEF/tencomp/tencomp.fdz")  # load a DataSet from file

# fd.viewer()  # start the viewer with no file opened
# fd.viewer(result)  # start the viewer and open the result DataSet
fd.viewer("simuEF/tencomp/tencomp.fdz")
