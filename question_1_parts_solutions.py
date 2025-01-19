import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from inspect import Traceback
from inspect import Traceback
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import  classification_report
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import classification_report,roc_auc_score,roc_curve
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score,precision_score,recall_score,classification_report
from ctypes import c_void_p
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import cross_val_score, KFold
from scipy.stats import ttest_ind




# Load the CSV file
file_path = "data.csv"  # Update with the actual path of your CSV file
data = pd.read_csv(file_path)

# Display the data
# print(data)

df=pd.DataFrame(data)
df['Tissue']=df['Tissue'].astype(str)
df.head()

le=LabelEncoder()
df['Tissues']=le.fit_transform(df['Tissue'])

for col in ['strand', 'CpG_Coordinates', '000', '001', '010', '011', '100', '101', '110', '111']:
    df[col] = df[col].apply(lambda x: ','.join(x) if isinstance(x, list) else str(x)).astype('category')

sns.countplot(x='000',hue='Tissues',data=df)
sns.countplot(x='001',hue='Tissues',data=df)
sns.countplot(x='010',hue='Tissues',data=df)
sns.countplot(x='011',hue='Tissues',data=df)
sns.countplot(x='100',hue='Tissues',data=df)
sns.countplot(x='101',hue='Tissues',data=df)
sns.countplot(x='110',hue='Tissues',data=df)
sns.countplot(x='111',hue='Tissues',data=df)
plt.title('PMP frequency by Tissue Type')
plt.show()

X=pd.get_dummies(df[['strand','CpG_Coordinates','000','001','010','011','100','101','110','111']],drop_first=True)
y=df['Tissues']
X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.2,random_state=0)

# model=LogisticRegression()
# model.fit(X_train,y_train)

# y_pred=model.predict(X_test)
# print(confusion_matrix(y_test,y_pred))
# print(classification_report(y_test,y_pred))
# print(accuracy_score(y_test,y_pred))


def calculate_coverage(df):
    binary_keys = ['000', '001', '010', '011', '100', '101', '110', '111']
    df['Coverage'] = df[binary_keys].sum(axis=1)
    return df

# Calculate coverage for the dataset
df = calculate_coverage(df)
# Display results
print(df[['strand', 'CpG_Coordinates', 'Coverage']].head())

coverage_data={
    'Tissues':['cfDNA','cfDNA','cfDNA','cfDNA','crDNA','crDNA','crDNA','crDNA'],
    'Cpg_Coverage':['1090','1090','1090','1090','1008','1003','1022','1009',]
}
df_coverage=pd.DataFrame(coverage_data)

# Convert 'Cpg_Coverage' to numeric before applying aggregation
df_coverage['Cpg_Coverage'] = pd.to_numeric(df_coverage['Cpg_Coverage'])

coverage_stats=df_coverage.groupby('Tissues')['Cpg_Coverage'].agg(
    Median='median',
    Mean='mean',
    Standard_Deviation='std',
)
coverage_stats['CV']=coverage_stats['Standard_Deviation']/coverage_stats['Mean']*100
print(coverage_stats)

sns.boxplot(x='Tissues',y='Cpg_Coverage',data=df_coverage)
plt.title('Distribution of CpG Coverage by Tissue Type')
plt.show()


sns.kdeplot(data=df_coverage,x='Cpg_Coverage',hue='Tissues',fill=True)
plt.title('Distribution of CpG Coverage by Tissue Type')
plt.show()

all_data = []
for entry in data:
    binary_keys = [key for key in entry if key.isdigit() and len(key) == 3]
    for i in range(len(entry['strand'])):
        # Compute coverage
        coverage = sum(int(entry[key][i]) for key in binary_keys)
        # Map tissue to 0 (cfDNA) or 1 (crDNA)
        tissue_label = 0 if entry['Tissue'][i] == 'cfDNA' else 1
        all_data.append({
            'Strand': entry['strand'][i],
            'CpG_Coordinates': entry['CpG_Coordinates'][i],
            'Coverage': coverage,
            'Tissue': entry['Tissue'][i],
            'Tissue_Label': tissue_label
        })

# Convert to DataFrame
df = pd.DataFrame(all_data)

# Group by tissue and compute statistics
statistics = df.groupby('Tissue')['Coverage'].agg(
    Median='median',
    Mean='mean',
    Standard_Deviation='std'
)
statistics['CV'] = statistics['Standard_Deviation'] / statistics['Mean'] * 100

# Output results
print(df)
print("\nStatistics:")
print(statistics)

methylated_states=['001','011','101','111']
unmethylated_states=['000','010','100','110']

all_data = []
for entery in data:
    binary_keys=[key for key in entery if key.isdigit() and len(key)==3]
    for i in range(len(entery['strand'])):
        total_count=sum(int(entery[key][i]) for key in binary_keys)
        methylated_count=sum(int(entery[key][i]) for key in binary_keys if key in methylated_states)
        pmp=(methylated_count/total_count)*100 if total_count>0 else 0
        all_data.append({
            'Strand':entery['strand'][i],
            'CpG_Coordinates':entery['CpG_Coordinates'][i],
            'Pmp':pmp,
            'Tissue': entery['Tissue'][i],  # Changed 'Tissues' to 'Tissue'
            'Total_count':total_count,
            'Methylated Count':methylated_count,
            'Tissue_Label':0 if entery['Tissue'][i]=='cfDNA' else 1
        })
df=pd.DataFrame(all_data)
pmp_by_tissue=df.groupby('Tissue')['Pmp'].mean().reset_index()

print("pmp data:")
print(df)
print("\npmp by tissue:")
print(pmp_by_tissue)

df['tissue_label']=df['Tissue'].map({'cfDNA':0,'crDNA':1})
X=df[['Pmp']]
y=df['tissue_label']
model=LogisticRegression()
model.fit(X,y)

df['Probability']=model.predict_proba(X)[:,1]
print(df)

df['P_value']=1-df['Probability']

print(classification_report(y,model.predict(X)))

roc_auc = roc_auc_score(y, model.predict_proba(X)[:, 1])
print(f'ROC AUC: {roc_auc}')


df['VRF']=df['Pmp']*df['P_value']
print(df)

df['Vrf']=df['Methylated Count']/df['Total_count']
mean_vrf= df.groupby('Tissue')['Vrf'].mean().reset_index()
print(mean_vrf)
# print(df)

vrf_cfdna=df[df['Tissue']=='cfDNA']['Vrf']
vrf_crdna=df[df['Tissue']=='crDNA']['Vrf']

t_statistic, p_value = ttest_ind(vrf_cfdna, vrf_crdna)

print(f"T-statistic: {t_statistic}")
print(f"P-value: {p_value}")

# Get the data from the original 'data' list
all_data_with_binary_cols = []
for entry in data:
    binary_keys = [key for key in entry if key.isdigit() and len(key) == 3]
    for i in range(len(entry['strand'])):
        # Get binary values for this row
        binary_values = [int(entry[key][i]) for key in binary_keys]

        # Append data for this row
        row_data = {
            'Strand': entry['strand'][i],
            'CpG_Coordinates': entry['CpG_Coordinates'][i],
            'Pmp': df.loc[(df['Strand'] == entry['strand'][i]) & (df['CpG_Coordinates'] == entry['CpG_Coordinates'][i]), 'Pmp'].iloc[0], # Get Pmp value from existing df
            'Tissue': entry['Tissue'][i],
            'Total_count': df.loc[(df['Strand'] == entry['strand'][i]) & (df['CpG_Coordinates'] == entry['CpG_Coordinates'][i]), 'Total_count'].iloc[0], # Get Total_count value from existing df
            'Methylated Count': df.loc[(df['Strand'] == entry['strand'][i]) & (df['CpG_Coordinates'] == entry['CpG_Coordinates'][i]), 'Methylated Count'].iloc[0], # Get Methylated Count value from existing df
            'Tissue_Label': 0 if entry['Tissue'][i] == 'cfDNA' else 1
        }
        row_data.update(dict(zip(binary_keys, binary_values))) # Add binary columns

        # Calculate and add Vrf for this row
        row_data['Vrf'] = row_data['Methylated Count'] / row_data['Total_count'] if row_data['Total_count'] > 0 else 0

        all_data_with_binary_cols.append(row_data)

# Create a new DataFrame with the binary columns
df_with_binary_cols = pd.DataFrame(all_data_with_binary_cols)


# Calculate mean and std coverage using the new DataFrame
df_with_binary_cols['mean_coverage'] = df_with_binary_cols[['000', '001', '010', '011', '100', '101', '110', '111']].mean(axis=1)
df_with_binary_cols['std_coverage'] = df_with_binary_cols[['000', '001', '010', '011', '100', '101', '110', '111']].std(axis=1)

# Print the results
print(df_with_binary_cols[['Pmp', 'Vrf', 'mean_coverage', 'std_coverage']].head())


# Instead of using 'df', use 'df_with_binary_cols' which contains the required columns.
X = df_with_binary_cols[['Vrf', 'mean_coverage', 'std_coverage']]
y = df_with_binary_cols['Tissue'].map({'cfDNA': 0, 'crDNA': 1})

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

rf_model = RandomForestClassifier(random_state=0)
rf_model.fit(X_train, y_train)
rf_predictions = rf_model.predict(X_test)

lr_model = LogisticRegression()
lr_model.fit(X_train, y_train)
lr_predictions = lr_model.predict(X_test)

print(classification_report(y_test, rf_predictions))
print(classification_report(y_test, lr_predictions))


importances = rf_model.feature_importances_
features= X.columns

plt.figure(figsize=(10, 6))
plt.bar(features,importances,color='skyblue')
plt.xlabel('Features')
plt.ylabel('Importance')
plt.title('Random Forest Feature Importance')
plt.xticks(rotation=45)
plt.show()

accuracy=accuracy_score(y_test,rf_predictions)
precision=precision_score(y_test,rf_predictions)
recall=recall_score(y_test,rf_predictions)

# specificity=sum((y_test==1)&(rf_predictions==1))/sum(rf_predictions==1)
# f1_score=2*(precision*recall)/(precision+recall)

specificity=sum((y_test==0)&(rf_predictions==0))/sum(rf_predictions==0)
f1_score=2*(precision*recall)/(precision+recall)

print(f"Accuracy: {accuracy}")
print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"Specificity: {specificity}")
print(f"F1 Score: {f1_score}")


# Call the function to get the coverage data
coverage_data = calculate_coverage(data)

# Convert coverage data to DataFrame
df = pd.DataFrame(coverage_data)

# Calculate median and CV for coverage by Tissue
coverage_stats = df.groupby('Tissue')['Coverage'].agg(
    Median='median',
    Mean='mean',
    Standard_Deviation='std'
)
coverage_stats['CV'] = coverage_stats['Standard_Deviation'] / coverage_stats['Mean'] * 100

# Display the results
print("\nCoverage Statistics:")
print(coverage_stats)

# Step 4: Plotting the coverage distribution by tissue
# Dynamically select all unique tissue types
tissue_types = df['Tissue'].unique()

# Prepare data for the boxplot
coverage_by_tissue = [df[df['Tissue'] == tissue]['Coverage'] for tissue in tissue_types]

# Create the boxplot
plt.figure(figsize=(10, 6))
plt.boxplot(coverage_by_tissue, labels=tissue_types)
plt.title('Coverage Distribution by Tissue Type')
plt.ylabel('Coverage')
plt.xticks(rotation=45)  # Rotate x-axis labels for better visibility
plt.show()

