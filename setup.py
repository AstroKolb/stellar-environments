import os

# open icons file
icons = open('icons')


# determine prefix
prefix = ''
for LINE in icons:

   if 'prefix' in LINE: prefix = LINE.split("'")[1]


# copy icons file
os.system('cp -p icons output/'+prefix+'.icons')
