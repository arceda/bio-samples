# Imports
from . import k_mers as K_mers
from . import matrices
from .  import classifiers

import matplotlib.pyplot as plt
from sklearn.metrics import f1_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold

def extractKmers(T, training_data, k_min, k_max, features_min, features_max, fig_file, kmers_file):
	# List of different lengths of k-mer
	k_mers_range = range(k_min, k_max + 1)
	# Features range
	features_range = range(features_min, features_max + 1)
	# Contains scores lists of different lengths of k
	scores_list = []
	# Contains list of k-mer for each iteration of rfe
	supports_list = []
	# Index to determine the number of k-mer
	index = 0
	# Best score of classification (F-measure)
	best_score = 0
	# Optimal score in relation with treshold
	optimal_score = 0
	# Best k-mer list
	best_k_mers = []
	# Best length of k
	best_k_length = 0
	# Classifier
	clf = classifiers.svm()
	# Splits number for evaluation
	n_splits = 5
	# Evaluation methode
	skf = StratifiedKFold(n_splits = n_splits, shuffle = False, random_state = None)
	# Corresponds to the percentage of features to remove at each iteration
	preliminary_rfe_step = 0.1

	# Perform the analysis for the different sizes of k
	for k in k_mers_range:

		print("\n\nBeginning of the " + str(k) + "_mer(s) analysis")

		# Generarte list of k-mer
		print("Generate K-mers...")
		k_mers = K_mers.generate_K_mers(training_data, k)
	
		# Genrate matrice attributes and matrice class
		print("Generate matrices...")
		X, y = matrices.generateXYMatrice(training_data, k_mers, k)

		# Identify maximum value of the matrice
		X_max = max(max(X))

		# Apply Min-max scaling or not
		if X_max > 1:
		
			print("Apply linearly scaling each attribute to the range [0, 1]")
			minMaxScaler = MinMaxScaler(feature_range=(0, 1), copy = False)
			X = minMaxScaler.fit_transform(X)
		else: print("Scaling not required ")
	

		# If more than features_max  apply RFE (remove 10 % of features to remove at each iteration)
		if len(X[0]) > features_max:
			print("Preliminary - RFE...")	
			rfe = RFE(estimator = clf, n_features_to_select = features_max, step = preliminary_rfe_step)
			X = rfe.fit_transform(X, y)

			# Update list of k_mers
			for i, value in enumerate(rfe.support_):
				if value == False: k_mers[i] = None
			k_mers = list(filter(lambda a: a != None, k_mers))


		# Classification score of each iteration
		scores = []
		# K-mer list used at each iteration
		supports = []

		# Encode label 
		labelEncodel = LabelEncoder()
		y = labelEncodel.fit_transform(y)

		for n in features_range:
			print("\rRFE :", round(n / features_max * 100, 0), "%", end='')
			f_measure = 0
			k_mers_rfe = []
			rfe = RFE(estimator = clf, n_features_to_select = n, step = 1)
			X_rfe = rfe.fit_transform(X, y)

			# Update list of k_mers
			for i, j in enumerate(rfe.support_):
				if j == True: k_mers_rfe.append(k_mers[i])
		
			# Evaluation of attributes with F-measure and Cross-validation
			for train_index, test_index in StratifiedKFold(n_splits = n_splits, shuffle=False, random_state=None).split(X_rfe, y):
				X_train, X_test = list(X_rfe[test_index]), list(X_rfe[train_index])
				y_train, y_test = list(y[test_index]), list(y[train_index])

				# Prediction
				clf.fit(X_train, y_train)
				y_pred = clf.predict(X_test)
			
				# Calcul metric scores
				f_measure = f_measure + f1_score(y_test, y_pred, average ="weighted")

			# Calcul mean F-measure
			mean_f_measure = f_measure / n_splits

			# Save scores
			scores.append(mean_f_measure)
			supports.append(k_mers_rfe)

		# Save list of score
		scores_list.append(scores)
		supports_list.append(supports)

	# Identify best solution
	for i, s in enumerate(scores_list):
		if max(s) > best_score:
			best_score = max(s)
			index = s.index(max(s))
			best_k_length = k_mers_range[i]
			best_k_mers = supports_list[i][index]
		elif max(s) == best_score:
			if s.index(max(s)) < index:
				best_score = max(s)
				index = s.index(max(s))
				best_k_length = k_mers_range[i]
				best_k_mers = supports_list[i][index]
		else: pass

	# Identify optimal solution
	for i, l in enumerate(scores_list):
		for j, s in enumerate(l):
			if s >=  best_score * T and j <= index: 
				optimal_score = s
				index = j
				best_k_length = k_mers_range[i]
				best_k_mers = supports_list[i][index]
				print("\nChange optimal solution")
			
	if optimal_score == 0: optimal_score = best_score

	# Save plot results
	fig = plt.figure(figsize=(12, 10))
	for i, s in enumerate(scores_list):
		label = str(k_mers_range[i]) + "-mers"
		plt.plot(features_range, s, label= label)
	plt.ylabel("F-measure")
	plt.xlabel("Number of features")
	plt.axvline(index + 1, linestyle=':', color='r')
	title = "F-measure : " + str(optimal_score) + " K-mer size : " + str(best_k_length) + " Number of features : " + str(index + 1)
	plt.title(title)
	plt.legend()
	print("Saving figure to " + fig_file + ".png")
	fname = str(fig_file)
	plt.savefig(fname)


	# Open a file
	print("Writing to " + kmers_file)
	f = open(kmers_file, "w")
	# Add k-mers
	for i in best_k_mers: f.write(str(i) + "\n");
	# Close opend file
	f.close()

	# Return identified k-mers and their length
	return best_k_mers, best_k_length
