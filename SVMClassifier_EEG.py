import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

#Load data
data_loc = 'C:/Users/varsh/Desktop/PhD/FinalPaper/'

data = pd.read_csv(data_loc + 'ThesisAnalysis_DataPresent_WithoutOutliers_5.29.csv')

#Creating age category -- Below and Above 40 years

data['age_group'] = data['Age'].apply(lambda age: 0 if age <= 40 else 1)

#%% Binding | SVM, Feature Extraction

binding_features = ['Coh12Amp4chan','Coh12Lat4chan',
                    'Coh20Amp4chan','Coh20Lat4chan',
                    'Incoh12P1Amp4chan','Incoh12P1Lat4chan',
                    'Incoh20P1Amp4chan','Incoh20P1Lat4chan',
                    'Incoh12P2Amp4chan','Incoh12P2Lat4chan',
                    'Incoh20P2Amp4chan','Incoh20P2Lat4chan',
                    'Coh12SS_4chan',
                    'Coh20SS_4chan',
                    'Incoh12SS_4chan',
                    'Incoh20SS_4chan']
                    # 'SSDiff12',
                    # 'SSDiff20']

# Select relevant columns
columns = ['age_group', 'Subject'] + binding_features
dat_binding = data[columns].dropna().reset_index(drop=True)

# Separate target (y) and features (X)
y_binding = dat_binding['age_group']
binding_model = dat_binding.drop(['age_group', 'Subject'], axis=1)

# Set parameters
n_runs = 10
n_features = binding_model.shape[1]
feature_names = binding_model.columns

# Initialize storage
binding_coef = np.zeros((n_runs, n_features))
binding_coef_abs = np.zeros((n_runs, n_features))
binding_confusion_matrices = []
binding_accuracy_scores = []
binding_cv_scores_all = []

# Define classifier
binding_clf = SVC(kernel='linear', C=1.0)

# Loop over multiple splits
for i in range(n_runs):

    # Train-test split
    binding_train, binding_test, y_train_binding, y_test_binding = train_test_split(
        binding_model, y_binding, test_size=0.1)

    # Fit classifier
    binding_clf.fit(binding_train, y_train_binding)

    # Store raw and absolute coefficients
    binding_coef[i] = binding_clf.coef_[0]
    binding_coef_abs[i] = np.abs(binding_clf.coef_[0])

    # Predict and evaluate
    pred_binding = binding_clf.predict(binding_test)
    binding_confusion_matrices.append(confusion_matrix(y_test_binding, pred_binding))
    binding_accuracy_scores.append(accuracy_score(y_test_binding, pred_binding))

    # Run cross-validation on the training split
    cv_scores = cross_val_score(binding_clf, binding_train, y_train_binding, cv=10)
    binding_cv_scores_all.append(cv_scores.mean())

# Compute averages
binding_coef_mean = binding_coef.mean(axis=0)
binding_coef_abs_mean = binding_coef_abs.mean(axis=0)
binding_confusion_matrix_mean = sum(binding_confusion_matrices) / len(binding_confusion_matrices)
# binding_confusion_matrix_mean = binding_confusion_matrix_mean.astype(int)
binding_confusion_matrix_mean = np.round(binding_confusion_matrix_mean).astype(int)
binding_accuracy_mean = np.mean(binding_accuracy_scores)
binding_cv_accuracy_mean = np.mean(binding_cv_scores_all)

# Print Summary Metrics ====
print("Mean Binding Accuracy Score:", binding_accuracy_mean)
print("Mean Binding CV Score:", binding_cv_accuracy_mean)
print("Mean Binding Confusion Matrix:", binding_confusion_matrix_mean)

#%% Plotting Binding SVM Results

# Feature Importance (Absolute Coefficients)
feature_names = binding_model.columns
sorted_idx = np.argsort(binding_coef_abs_mean)[::-1]

plt.figure(figsize=(8, 6))
plt.bar(range(n_features), binding_coef_abs_mean[sorted_idx], color='lightsteelblue')
plt.xticks(range(n_features), feature_names[sorted_idx], rotation=45, ha='right')
plt.ylabel("Mean Absolute Coefficient", fontsize = 16, fontweight = "bold")
plt.title("Binding Features (Linear SVM)", fontsize = 18, fontweight = "bold")

# Add text for mean accuracy and CV score on the plot
mean_accuracy_text = f"Mean Accuracy: {binding_accuracy_mean:.3f}"
mean_cv_score_text = f"Mean CV Score: {binding_cv_accuracy_mean:.3f}"

plt.text(1.5, max(binding_coef_abs_mean) * 0.9, mean_accuracy_text, fontsize=14, color='black')
plt.text(1.5, max(binding_coef_abs_mean) * 0.85, mean_cv_score_text, fontsize=14, color='black')

plt.tight_layout()
plt.show()
plt.savefig(data_loc + 'ACFeature_FeatureImportance_4.18.png', dpi=500)

# Feature Directionality (Raw Coefficients)
plt.figure(figsize=(8, 6))
plt.bar(range(n_features), binding_coef_mean[sorted_idx], color='slateblue')
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.xticks(range(n_features), feature_names[sorted_idx], rotation=45, ha='right')
plt.ylabel("Mean Coefficient", fontsize = 16, fontweight = "bold")
plt.title("Binding Features – Directional Weights (Linear SVM)", fontsize = 18, fontweight = "bold")
plt.tight_layout()
plt.show()

# Mean Confusion Matrix
plt.figure(figsize=(5, 4))
sns.heatmap(binding_confusion_matrix_mean, annot=True, cmap='Blues')
plt.title("Mean Confusion Matrix (Binding Model)")
plt.xlabel("Predicted Age Group")
plt.ylabel("True Age Group")
plt.tight_layout()
plt.show()

# Boxplot of Feature Distributions by Age Group
dat_binding_melt = dat_binding.melt(id_vars='age_group', value_vars=binding_features)

plt.figure(figsize=(9, 6))
sns.boxplot(x='variable', y='value', hue='age_group', data=dat_binding_melt)
handles, labels = plt.gca().get_legend_handles_labels()
labels = ['<=40' if label == '0' else '>40' for label in labels]
plt.legend(handles, labels, title='Age Group')
plt.xticks(rotation=45)
plt.xlabel ('')
plt.ylabel("Coefficient", fontsize = 16, fontweight = "bold")
plt.title('Distribution of Binding Features by Age Group', fontsize = 18, fontweight = "bold")
plt.tight_layout()
plt.show()

plt.savefig(data_loc + 'ACFeature_FeatureImp_AgeGroup_4.18.png', dpi=500)

# Save Single Weighted Feature for each subject
# Weighted sum of features using abs and not abs mean coefficients
abs_binding = binding_coef_abs_mean * binding_model[feature_names]
abs_binding_SF = abs_binding.sum(axis=1)

binding = binding_coef_mean * binding_model[feature_names]
binding_SF = binding.sum(axis=1)

binding_all = pd.concat([dat_binding[['Subject']], abs_binding_SF.rename('AC_Feature_Abs'), binding_SF.rename('AC_Feature_NotAbs')], axis=1)

binding_all.to_csv(data_loc + 'Binding_SF_Evoked_6.10.csv', index = False)

#%% GDT | SVM

# Setting up...
GDT_features = ['GapP14chanAmp16','GapP14chanLat16',
                'GapP24chanAmp16','GapP24chanLat16',
                'GapP14chanAmp32','GapP14chanLat32',
                'GapP24chanAmp32','GapP24chanLat32',
                'GapP14chanAmp64','GapP14chanLat64',
                'GapP24chanAmp64','GapP24chanLat64',
                'SS16', 'SS32', 'SS64']

GDT_columns =  ['age_group', 'Subject'] + GDT_features
dat_GDT = data[GDT_columns].dropna().reset_index(drop=True)

# Separate target (y) and features (X)
y_GDT = dat_GDT['age_group']
GDT_model = dat_GDT.drop(['age_group', 'Subject'], axis=1)

feature_names_GDT = GDT_model.columns
n_features_GDT = GDT_model.shape[1]
n_runs = 10

# Initialize storage
GDT_coef = np.zeros((n_runs, n_features_GDT))
GDT_coef_abs = np.zeros((n_runs, n_features_GDT))
GDT_confusion_matrices = []
GDT_accuracy_scores = []
GDT_cv_scores_all = []

# Classifier
GDT_clf = SVC(kernel='linear', C=1.0, random_state=2)

# Run model n times
for i in range(n_runs):

    # Train-test split
    GDT_train, GDT_test, y_train_GDT, y_test_GDT = train_test_split(
        GDT_model, y_GDT, test_size=0.1)

    # Fit classifier
    GDT_clf.fit(GDT_train, y_train_GDT)

    # Store coefficients
    GDT_coef[i] = GDT_clf.coef_[0]
    GDT_coef_abs[i] = np.abs(GDT_clf.coef_[0])

    # Test performance
    pred_GDT = GDT_clf.predict(GDT_test)
    GDT_confusion_matrices.append(confusion_matrix(y_test_GDT, pred_GDT))
    GDT_accuracy_scores.append(accuracy_score(y_test_GDT, pred_GDT))

    # Cross-validation score on training data
    cv_scores = cross_val_score(GDT_clf, GDT_train, y_train_GDT, cv=8)
    GDT_cv_scores_all.append(cv_scores.mean())

# Summary metrics
GDT_coef_mean = GDT_coef.mean(axis=0)
GDT_coef_abs_mean = GDT_coef_abs.mean(axis=0)
GDT_confusion_matrix_mean = sum(GDT_confusion_matrices) / len(GDT_confusion_matrices)
# GDT_confusion_matrix_mean = GDT_confusion_matrix_mean.astype(int)
GDT_confusion_matrix_mean = np.round(GDT_confusion_matrix_mean).astype(int)
GDT_accuracy_mean = np.mean(GDT_accuracy_scores)
GDT_cv_accuracy_mean = np.mean(GDT_cv_scores_all)

print("Mean GDT Accuracy Score:", GDT_accuracy_mean)
print("Mean GDT CV Score:", GDT_cv_accuracy_mean)
print("Mean GDT Confusion Matrix:", GDT_confusion_matrix_mean)

#%% Plotting GDT Results

# Feature Importance (Absolute Coefficients)
feature_names_GDT = GDT_model.columns
sorted_idx = np.argsort(GDT_coef_abs_mean)[::-1]

plt.figure(figsize=(8, 6))
plt.bar(range(len(feature_names_GDT)), GDT_coef_abs_mean[sorted_idx], color='mediumseagreen')
plt.xticks(range(len(feature_names_GDT)), feature_names_GDT[sorted_idx], rotation=45, ha='right')
plt.ylabel("Mean Absolute Coefficient", fontsize = 16, fontweight = "bold")
plt.title("GDT Features (Linear SVM)", fontsize = 18, fontweight = "bold")

# Add text for mean accuracy and CV score on the plot
mean_accuracy_text = f"Mean Accuracy: {GDT_accuracy_mean:.3f}"
mean_cv_score_text = f"Mean CV Score: {GDT_cv_accuracy_mean:.3f}"

plt.text(1, max(GDT_coef_abs_mean) * 0.9, mean_accuracy_text, fontsize=14, color='black')
plt.text(1, max(GDT_coef_abs_mean) * 0.85, mean_cv_score_text, fontsize=14, color='black')

plt.tight_layout()
plt.show()

plt.savefig(data_loc + 'WCFeature_FeatureImportance_4.18.png', dpi=500)

# Feature Directionality (Raw Coefficients)
# plt.figure(figsize=(10, 6))
# plt.bar(range(n_features), binding_coef_mean[sorted_idx], color='slateblue')
# plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
# plt.xticks(range(n_features), feature_names[sorted_idx], rotation=45, ha='right')
# plt.ylabel("Mean Coefficient")
# plt.title("Binding Features – Directional Weights (Linear SVM)")
# plt.tight_layout()
# plt.show()

# Mean Confusion Matrix
plt.figure(figsize=(5, 4))
sns.heatmap(GDT_confusion_matrix_mean, annot=True, cmap='Blues')
plt.title("Mean Confusion Matrix (GDT Model)")
plt.xlabel("Predicted Age Group")
plt.ylabel("True Age Group")
plt.tight_layout()
plt.show()

# Boxplot of Feature Distributions by Age Group
dat_GDT_melt = dat_GDT.melt(id_vars='age_group', value_vars=GDT_features)

plt.figure(figsize=(9, 6))
sns.boxplot(x='variable', y='value', hue='age_group', data=dat_GDT_melt)

handles, labels = plt.gca().get_legend_handles_labels()
labels = ['<=40' if label == '0' else '>40' for label in labels]
plt.legend(handles, labels, title='Age Group')
plt.xticks(rotation=45)
plt.xlabel ('')
plt.ylabel("Coefficient", fontsize = 16, fontweight = "bold")
plt.title('Distribution of GDT Features by Age Group', fontsize = 18, fontweight = "bold")
plt.tight_layout()
plt.show()

plt.savefig(data_loc + 'WCFeature_FeatureImp_AgeGroup_4.18.png', dpi=500)

# Save Single Weighted Feature Per Subject
# Weighted sum of features using abs mean coefficients
abs_GDT = GDT_coef_abs_mean * GDT_model[feature_names_GDT]
abs_GDT_SF = abs_GDT.sum(axis=1)

GDT = GDT_coef_mean * GDT_model[feature_names_GDT]
GDT_SF = GDT.sum(axis=1)

GDT_all = pd.concat([dat_GDT[['Subject']], abs_GDT_SF.rename('WC_Feature_Abs'), GDT_SF.rename('WC_Feature_NotAbs')], axis=1)

GDT_all.to_csv(data_loc + 'GDT_SF_Evoked_6.10.csv', index = False)

#%% Comparing the models

# 1. Plotting Comparison of Accuracy and CV Scores for Binding SVM and GDT SVM

# Compute standard deviations
binding_accuracy_std = np.std(binding_accuracy_scores)
GDT_accuracy_std = np.std(GDT_accuracy_scores)

binding_cv_std = np.std(binding_cv_scores_all)
GDT_cv_std = np.std(GDT_cv_scores_all)

# Plotting
plt.figure(figsize=(8, 6))

models = ['Binding SVM', 'GDT SVM']
accuracy_scores = [binding_accuracy_mean, GDT_accuracy_mean]
cv_scores = [binding_cv_accuracy_mean, GDT_cv_accuracy_mean]
accuracy_std = [binding_accuracy_std, GDT_accuracy_std]
cv_std = [binding_cv_std, GDT_cv_std]

bar_width = 0.35
index = np.arange(len(models))

# Accuracy bars with SD
plt.bar(index, accuracy_scores, bar_width, yerr=accuracy_std,
        label='Accuracy', color='lightblue', capsize=5)

# # CV Score bars with SD
plt.bar(index + bar_width, cv_scores, bar_width, yerr=cv_std,
        label='CV Score', color='salmon', capsize=5)

# Add text for mean accuracy and CV score on the plot
mean_accuracy_binding = f"Binding Mean Accuracy: {binding_accuracy_mean:.3f}"
mean_cv_binding = f"Binding Mean CV Score: {binding_cv_accuracy_mean:.3f}"
mean_accuracy_GDT = f"GDT Mean Accuracy: {GDT_accuracy_mean:.3f}"
mean_cv_GDT = f"GDT Mean CV Score: {GDT_cv_accuracy_mean:.3f}"

plt.text(0, 0.95, mean_accuracy_binding, fontsize=12, color='black')
plt.text(0, 0.9, mean_cv_binding, fontsize=12, color='black')
plt.text(1, 0.95, mean_accuracy_GDT, fontsize=12, color='black')
plt.text(1, 0.9, mean_cv_GDT, fontsize=12, color='black')


plt.ylabel('Score', fontsize=12)
plt.xticks(index + bar_width / 2, models, rotation=0)

plt.ylim(0, 1.1)
plt.legend()
plt.tight_layout()
plt.show()

plt.savefig(data_loc + 'Binding_vs_GDT_Accuracy,CVScore_Comparison_4.18.png', dpi=500)

#%% Comparing the coefficients of both the models

plt.figure(figsize=(10, 6))

binding_sorted_idx = np.argsort(binding_coef_abs_mean)[::-1]  # Sort coefficients for Binding SVM
GDT_sorted_idx = np.argsort(GDT_coef_abs_mean)[::-1]  # Sort coefficients for GDT SVM

# Get the top N features to compare (for visualization purposes)
top_n = 10
top_binding_coef = binding_coef_abs_mean[binding_sorted_idx[:top_n]]
top_GDT_coef = GDT_coef_abs_mean[GDT_sorted_idx[:top_n]]

top_binding_features = np.array(feature_names)[binding_sorted_idx[:top_n]]
top_GDT_features = np.array(feature_names_GDT)[GDT_sorted_idx[:top_n]]

bar_width = 0.5

plt.bar(np.arange(top_n), top_binding_coef, bar_width, label='Binding SVM', color='lightsteelblue')

plt.bar(np.arange(top_n, top_n + top_n), top_GDT_coef, bar_width, label='GDT SVM', color='mediumseagreen')

plt.xlabel('Features', fontsize=12)
plt.ylabel('Mean Absolute Coefficient', fontsize=12)
plt.title('Comparison of First 10 Coefficients | SVM', fontsize=14)

combined_features = np.concatenate((top_binding_features, top_GDT_features))

plt.xticks(np.arange(0, 2 * top_n), combined_features, rotation=45, ha='right')

mean_accuracy_binding = f"Binding Mean Accuracy: {binding_accuracy_mean:.3f}"
mean_cv_binding = f"Binding Mean CV Score: {binding_cv_accuracy_mean:.3f}"
mean_accuracy_GDT = f"GDT Mean Accuracy: {GDT_accuracy_mean:.3f}"
mean_cv_GDT = f"GDT Mean CV Score: {GDT_cv_accuracy_mean:.3f}"

plt.text(1, max(binding_coef_abs_mean) * 0.9, mean_accuracy_binding, fontsize=12, color='darkblue')
plt.text(1, max(binding_coef_abs_mean) * 0.85, mean_cv_binding, fontsize=12, color='darkblue')
plt.text(10, max(binding_coef_abs_mean) * 0.9, mean_accuracy_GDT, fontsize=12, color='darkgreen')
plt.text(10, max(binding_coef_abs_mean) * 0.85, mean_cv_GDT, fontsize=12, color='darkgreen')

plt.legend()
plt.tight_layout()
plt.show()

plt.savefig(data_loc + 'Binding_vs_GDT_Coeff_Comparison_4.18.png', dpi=500)