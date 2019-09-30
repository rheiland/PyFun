from pyMCDS import pyMCDS
import numpy as np

#mcds = pyMCDS('output00000600.xml', '/Users/heiland/Documents/Paul/John/leader_follower_invasionSep2019/frozen')
#mcds = pyMCDS('output00000600.xml', '.')
fname = 'output00000600.xml'
print('reading ',fname)
mcds1 = pyMCDS(fname)

print(mcds1.get_time())
# print(mcds1.get_menv_species_list())
#print(mcds1.get_concentrations('quorum'))

cell_ids = mcds1.data['discrete_cells']['ID']
print(cell_ids.shape)
print(cell_ids[:4])

#cell_vec = np.zeros((cell_ids.shape, 3))
cell_vec = np.zeros((cell_ids.shape[0], 3))
vec_list = ['position_x', 'position_y', 'position_z']
for i, lab in enumerate(vec_list):
  cell_vec[:, i] = mcds1.data['discrete_cells'][lab]
xvals = cell_vec[:, 0]
print('x range: ',xvals.min(), xvals.max())
yvals = cell_vec[:, 1]
print('y range: ',yvals.min(), yvals.max())

dist2 = xvals**2 + yvals**2
dist = dist2**0.5
print("avg dist of all cells = ",dist.sum()/len(dist))

#print(mcds1.get_cell_variables())
