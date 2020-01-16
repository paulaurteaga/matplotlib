# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import sem
import numpy as np
from scipy.stats import linregress

# Study data files
mouse_metadata = "data/Mouse_metadata.csv"
study_results = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata)
study_results = pd.read_csv(study_results)

# Combine the data into a single dataset
merged=pd.merge(mouse_metadata,study_results, on='Mouse ID')
print(merged.head())
# Generate a summary statistics table of mean, median, ariance, standard deviation, and SEM of the tumor volume for each regimen
grouped=merged.groupby('Drug Regimen')
mean=grouped['Tumor Volume (mm3)'].mean()
median=grouped['Tumor Volume (mm3)'].median()
variance=grouped['Tumor Volume (mm3)'].var()
stDev=grouped['Tumor Volume (mm3)'].std()
sem=grouped['Tumor Volume (mm3)'].sem()
dic={'Mean': mean, 'median': median, 'variance': variance, 'Standar devition': stDev, 'SEM': sem}
summary=pd.DataFrame(dic)
print(summary)
# Generate a bar plot showing number of data points for each treatment regimen using pandas
number_points=grouped['Drug Regimen'].count()
regimens=summary.index
dic={'Regimens': regimens, "Data points": number_points}
graph=pd.DataFrame(dic)
ax=graph.plot(kind='bar', figsize=(10,7), rot=0,width=0.75,title='Number of data points per treatment regimen', legend=True)
ax.set_xlabel("Treatment regimens")
ax.set_ylabel("Number of data points")
plt.show(ax)
# Generate a bar plot showing number of data points for each treatment regimen using pyplot
number_points=grouped['Drug Regimen'].count()
regimens=summary.index
plt.figure(figsize=(10,7))
plt.bar(regimens, number_points)
plt.title('Number of data points per treatment regimen')
plt.xlabel("Treatment regimens")
plt.ylabel("Number of data points")
plt.legend(['Data points'])
plt.show()
# Generate a pie plot showing the distribution of female versus male mice using pandas
grouped=merged.groupby('Sex')
count=grouped['Mouse ID'].nunique()
ax=count.plot(kind='pie', legend=True, title='Female vs Male mice distribution', autopct='%1.1f%%')
ax.legend(loc='lower right')
plt.show(ax)
# Generate a pie plot showing the distribution of female versus male mice using pyplot
plt.pie(count,labels=count.index, autopct='%1.1f%%')
plt.title('Female vs Male mice distribution')
plt.legend(count.index, loc='lower right')
plt.show()


# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
max_tumor = merged.groupby(["Mouse ID"]).max()
max_tumor = max_tumor.reset_index()
merged_data = max_tumor[['Mouse ID','Timepoint']].merge(merged,on=['Mouse ID','Timepoint'],how="left")
capomulin = merged_data.loc[merged_data["Drug Regimen"] == "Capomulin"]['Tumor Volume (mm3)'] 
ramicane = merged_data.loc[merged_data["Drug Regimen"] == "Ramicane"]['Tumor Volume (mm3)'] 
infubinol = merged_data.loc[merged_data["Drug Regimen"] == "Infubinol"]['Tumor Volume (mm3)'] 
ceftamin = merged_data.loc[merged_data["Drug Regimen"] == "Ceftamin"]['Tumor Volume (mm3)'] 
#change column orders
max_tumor=max_tumor[['Mouse ID','Timepoint', 'Tumor Volume (mm3)', 'Metastatic Sites','Drug Regimen', 'Sex', 'Age_months', 'Weight (g)']]
# Calculate the IQR and quantitatively determine if there are any potential outliers. grouped=merged.groupby(['Drug Regimen','Mouse ID'])
#First Capomulin regimen
cap_quartiles = capomulin.quantile([.25,.5,.75])
cap_lowerq = cap_quartiles[0.25]
cap_upperq = cap_quartiles[0.75]
cap_iqr = cap_upperq-cap_lowerq
cap_lower_bound = cap_lowerq - (1.5*cap_iqr)
cap_upper_bound = cap_upperq + (1.5*cap_iqr)
print(f"The lower quartile of avg tumor volume is: {round(cap_lowerq,2)}")
print(f"The upper quartile of avg tumor volume is: {round(cap_upperq,2)}")
print(f"The interquartile range of avg tumor volume is: {round(cap_iqr,2)}")
print(f"The the median of avg tumor volume is: {round(cap_quartiles[0.5],2)} ")

lower_bound = cap_lowerq - (1.5*cap_iqr)
upper_bound = cap_upperq + (1.5*cap_iqr)
print(f"Values below {round(cap_lower_bound,2)} could be outliers.")
print(f"Values above {round(cap_upper_bound,2)} could be outliers.")
# Calculate the IQR and quantitatively determine if there are any potential outliers
#Second ramicane regimen
ram_quartiles = ramicane.quantile([.25,.5,.75])
ram_lowerq = ram_quartiles[0.25]
ram_upperq = ram_quartiles[0.75]
ram_iqr = ram_upperq-ram_lowerq
ram_lower_bound = ram_lowerq - (1.5*ram_iqr)
ram_upper_bound = ram_upperq + (1.5*ram_iqr)
print(f"The lower quartile of avg tumor volume is: {round(ram_lowerq,2)}")
print(f"The upper quartile of avg tumor volume is: {round(ram_upperq,2)}")
print(f"The interquartile range of avg tumor volume is: {round(ram_iqr,2)}")
print(f"The the median of avg tumor volume is: {round(ram_quartiles[0.5],2)} ")

lower_bound = ram_lowerq - (1.5*ram_iqr)
upper_bound = ram_upperq + (1.5*ram_iqr)
print(f"Values below {round(ram_lower_bound,2)} could be outliers.")
print(f"Values above {round(ram_upper_bound,2)} could be outliers.")
# Calculate the IQR and quantitatively determine if there are any potential outliers. 
#third infubinol regimen
inf_quartiles = infubinol.quantile([.25,.5,.75])
inf_lowerq = inf_quartiles[0.25]
inf_upperq = inf_quartiles[0.75]
inf_iqr = inf_upperq-inf_lowerq
inf_lower_bound = inf_lowerq - (1.5*inf_iqr)
inf_upper_bound = inf_upperq + (1.5*inf_iqr)
print(f"The lower quartile of avg tumor volume is: {round(inf_lowerq,2)}")
print(f"The upper quartile of avg tumor volume is: {round(inf_upperq,2)}")
print(f"The interquartile range of avg tumor volume is: {round(inf_iqr,2)}")
print(f"The the median of avg tumor volume is: {round(inf_quartiles[0.5],2)} ")

lower_bound = inf_lowerq - (1.5*inf_iqr)
upper_bound = inf_upperq + (1.5*inf_iqr)
print(f"Values below {round(inf_lower_bound,2)} could be outliers.")
print(f"Values above {round(inf_upper_bound,2)} could be outliers.")
# Calculate the IQR and quantitatively determine if there are any potential outliers.
#fourth ceftamin regimen
cef_quartiles = ceftamin.quantile([.25,.5,.75])
cef_lowerq = cef_quartiles[0.25]
cef_upperq = cef_quartiles[0.75]
cef_iqr = cef_upperq-cef_lowerq
cef_lower_bound = cef_lowerq - (1.5*cef_iqr)
cef_upper_bound = cef_upperq + (1.5*cef_iqr)
print(f"The lower quartile of avg tumor volume is: {round(cef_lowerq,2)}")
print(f"The upper quartile of avg tumor volume is: {round(cef_upperq,2)}")
print(f"The interquartile range of avg tumor volume is: {round(cef_iqr,2)}")
print(f"The the median of avg tumor volume is: {round(cef_quartiles[0.5],2)} ")

lower_bound = cef_lowerq - (1.5*cef_iqr)
upper_bound = cef_upperq + (1.5*cef_iqr)
print(f"Values below {round(cef_lower_bound,2)} could be outliers.")
print(f"Values above {round(cef_upper_bound,2)} could be outliers.")
# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
outlier= dict(markerfacecolor='red',markersize=12)
plt.boxplot([capomulin,ramicane,infubinol,ceftamin],labels=['Capomulin','Ramicane','Infubinol','Ceftamin'],flierprops=outlier)
plt.ylabel('Final Tumor Volume (mm3)')
plt.show()
#Getting the possible Infubonil outlier
inf_outlier=0
for x in infubinol.iteritems():
    if x[1]<inf_lower_bound:
      inf_outlier=(x[1])  
print(f"The value {round(inf_outlier,2)} could be an outlier.")
# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
#First we just take the rows where the treatment is Capomulin by using .loc
capomulin_frame=merged.loc[merged['Drug Regimen']== 'Capomulin']
#We create a Series containing all the unique mouse ids for Capomulin treatments
mouse_ids=capomulin_frame['Mouse ID'].unique()
print(mouse_ids)
#We ask the user to select one of the mouses treated with Capomulin
selection=input('Select one of the mouses ID from above ')
#And we just get the timepoints and tumor volume info of the mouse id selected by the user
time_points=capomulin_frame.loc[capomulin_frame['Mouse ID']== selection, ['Timepoint']]
tumor_volume=capomulin_frame.loc[capomulin_frame['Mouse ID']== selection, ['Tumor Volume (mm3)']]
plt.figure(figsize=(8,5))
plt.plot(time_points,tumor_volume)
plt.title(f'Time points vs tumor volume of mouse ID number {selection}')
plt.xlabel('Time points')
plt.ylabel('Tummor volume in mm3')
plt.show()
# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
capomulin_frame=merged.loc[merged['Drug Regimen']== 'Capomulin']
capomulin_grouped=capomulin_frame.groupby('Mouse ID').mean()
weight=capomulin_grouped['Weight (g)']
volume=capomulin_grouped['Tumor Volume (mm3)']
plt.figure(figsize=(8,5))
plt.scatter(weight,volume,marker="+", facecolors="blue")
plt.title('Avg mouse weight vs Avg tumor volume')
plt.xlabel('Weight in grams')
plt.ylabel('Tummor volume in mm3')
plt.show()
# Calculate the correlation coefficient and linear regression model for mouse weight and average tumor volume for the Capomulin regimen
correlation = st.pearsonr(weight,volume)
print(f"The correlation between both factors is {round(correlation[0],2)}")
(slope, intercept, rvalue, pvalue, stderr) = linregress(weight,volume)
regress_values = weight * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.figure(figsize=(8,5))
plt.scatter(weight,volume,marker="+", color='blue')
plt.title('Avg mouse weight vs Avg tumor volume')
plt.xlabel('Weight in grams')
plt.ylabel('Tummor volume in mm3')
plt.plot(weight,regress_values,"r-")
plt.annotate(line_eq,(20,10),fontsize=15,color="red")
plt.show()
print(f'The ecuation is y={str(round(slope,2))}x + {str(round(intercept,2))}')