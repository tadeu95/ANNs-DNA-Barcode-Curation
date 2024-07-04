# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from keras import Sequential
from keras.layers import Dense
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from keras.callbacks import History
from keras.layers import Dropout
import joblib

#np.random.seed(42)
#tf.random.set_seed(42)

#Replace the file path with the actual path to your training set (tsv format)
file_path = 'C:/Users/username/files/labelled_for_nn_ready.tsv'
training_df = pd.read_csv(file_path, sep='\t')

#Extract features and target variable
counting_features = ['frequency_species_nn', 'frequency_bin_nn', 'num_this_bin_nn', 'unique_identifiers_nn', 'unique_institutions_nn','species_per_bin_nn','bin_per_species_nn']
numerical_features = ['percentage_sequences_in_cluster_nn']
categorical_features = ['species_nn', 'BIN_nn', 'genus_nn', 'ambiguous_name_nn', 'ingroup_synonym_nn']
target_variable = 'label_binary_nn'

X_counting = training_df[counting_features]
X_numerical = training_df[numerical_features]
X_categorical = training_df[categorical_features]
y = training_df[target_variable]

X_combined = pd.concat([X_counting, X_numerical, X_categorical], axis=1)

#Standardize the percentage feature
scaler = StandardScaler()
X_numerical_scaled = pd.DataFrame(scaler.fit_transform(X_numerical), columns=numerical_features)
scaler_path = 'C:/Users/username/files/scaler.save'
joblib.dump(scaler, scaler_path)

#Combine the features with standardized numerical feature
X_combined_scaled = pd.concat([X_numerical_scaled, np.log1p(X_counting), X_categorical], axis=1)

#Split the training data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(X_combined_scaled, y, test_size=0.1, random_state=32)

#Define the model
model = Sequential()
model.add(Dense(64, input_dim=X_combined_scaled.shape[1], activation='relu'))  
model.add(Dropout(0.05,seed=32))
model.add(Dense(32, activation='relu'))  
model.add(Dense(1, activation='sigmoid'))
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
history_callback = History()

#Train the model on the training set and validate using the validation set
history = model.fit(X_train, y_train, epochs=50, batch_size=64, validation_data=(X_val, y_val), callbacks=[history_callback])

#Extract featuress from the testing set
testing_file_path = 'C:/Users/username/files/testing_set_for_nn.tsv'
testing_df = pd.read_csv(testing_file_path, sep='\t')

X_test_counting = testing_df[counting_features]
X_test_numerical = testing_df[numerical_features]
X_test_categorical = testing_df[categorical_features]

#Combine the features with standardized numerical features and counting features
X_test_combined_scaled = pd.concat([pd.DataFrame(scaler.transform(X_test_numerical), columns=numerical_features), np.log1p(X_test_counting), X_test_categorical], axis=1)

#Make predictions on the testing set
predictions = model.predict(X_test_combined_scaled)

binary_predictions = (predictions > 0.5).astype(int)

#Load the ground truth labels from the testing set
ground_truth_labels = testing_df['ground_truth_nn']

#Evaluate the model's performance using the ground truth labels
accuracy = accuracy_score(ground_truth_labels, binary_predictions)
conf_matrix = confusion_matrix(ground_truth_labels, binary_predictions)
class_report = classification_report(ground_truth_labels, binary_predictions, digits=4)  # You can change 4 to any number of decimal places you

#Print evaluation metrics
print(f"Accuracy: {accuracy}")
print(f"Confusion Matrix:\n{conf_matrix}")
print(f"Classification Report:\n{class_report}")

#Create a dataframe with the original features and the predicted labels
result_df = pd.concat([testing_df, pd.Series(binary_predictions.flatten(), name='predicted_label_binary_nn')], axis=1)

#Save evaluation metrics
metrics_dict = {
    'Accuracy': accuracy,
    'Confusion Matrix': conf_matrix,
    'Classification Report': class_report
}
metrics_df = pd.DataFrame(metrics_dict.items(), columns=['Metric', 'Value'])

#Get the weights of the first layer
layer1_weights = model.layers[0].get_weights()[0]

#Get the feature names
feature_names = X_combined_scaled.columns

#Print the weights for each feature
for i, feature in enumerate(feature_names):
    print(f"Feature: {feature}")
    print(f"Weights: {layer1_weights[:, i]}\n")
    
#Create a confusion matrix
conf_matrix = confusion_matrix(ground_truth_labels, binary_predictions)

print(model.history.history.keys())

#Plot training and validation accuracy values
plt.plot(history.history['accuracy'], label='Training Accuracy')
plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.legend()
plt.title('Model Accuracy')
plt.savefig('C:/Users/username/files/training_validation_accuracy.png', dpi=300)
plt.show()

#Plot training and validation loss values
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.title('Model Loss')
plt.savefig('C:/Users/username/files/training_validation_loss.png', dpi=300)
plt.show()

#Save the result file with the model's predictions
result_file_name = 'C:/Users/username/files/result_file.tsv'
result_df.to_csv(result_file_name, sep='\t', index=False)

#Save the evaluation metrics in a csv file
metrics_file_name = 'C:/Users/username/files/evaluation_metrics.csv'
metrics_df.to_csv(metrics_file_name, index=False)

#Save the model in a h5 file 
model_file_name = 'C:/Users/username/files/my_model.h5'
model.save(model_file_name)

#Save the history of the model in a pkl file
history_file_name = 'C:/Users/username/files/history.pkl'
joblib.dump(history.history, history_file_name)

