
## Python packages used: Pridepy ppx,matplotlib,numpy

import pprint
from pridepy import Project
project = Project()
results = project.search_by_keywords_and_filters(
"taxonomy:Bordetella pertussis", # search for a specific species
"filter", # you can define a filter here, but you don't need to
1000, # maximum number of results
0, # pages - no need to change this
10, # date gap - no need to change this
"DESC", # order in which results are sorted
"submission_date", # sorting criterium
)
instruments = []

identifier = set()
for n, dataset in enumerate(results["_embedded"]["compactprojects"]):
 for i in dataset['instruments']:
    instruments.append(i)
   # instruments.append (str (dataset['instruments']))
    identifier.add(dataset["accession"])
    pprint.pprint(dataset["accession"])
    pprint.pprint(dataset["instruments"])  #Prints instruments used
    
    #print(f"{val}")
print(n)  # Prints the total number of datasets
pprint.pprint(identifier)

## Creating dictionary for counting instrunets for each dataset
Dict = {}
for i in instruments:
    val= i
    if(val in Dict):
          Dict[val]+= 1
    else:
        Dict[val]= 1
print(Dict)
    
## Creating histogram for the different instruments used        
from matplotlib import pyplot as plt      
  
  ## dataset
data = [instruments,val]
 
 # Plotting
plt.hist(data,bins=15,edgecolor = 'black')
plt.xlabel('Instruments')
plt.ylabel('Counts')
plt.title('Histogram of Instruments')
# plt.legend(instruments)
plt.show()   ## Shows plot
  




