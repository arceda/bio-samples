from . import matrices
from . import classifiers

from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict

def cross_validation(training_data, best_k_mers):
 
    # Generate  matrices
    best_k_length = len(best_k_mers[0])
    X, y = matrices.generateXYMatrice(training_data, best_k_mers, best_k_length)

	# Realize evaluation with CV + Classifier + Metrics
    clf = classifiers.svm()
    skf = StratifiedKFold(n_splits = 10, shuffle = False, random_state = None)
    y_pred = cross_val_predict(clf, X, y, cv = skf, n_jobs = 4)

    # Print results
    print("\nClassification report of model evaluation\n")
    print(classification_report(y, y_pred, digits = 3))
    print("Accuracy :", accuracy_score(y, y_pred) * 100, "%\n")
    print("Confusion matrix \n", confusion_matrix(y, y_pred))


def train_model(training_data, best_k_mers):

    # Generate  matrices
    best_k_length = len(best_k_mers[0])
    X_train, y_train = matrices.generateXYMatrice(training_data, best_k_mers, best_k_length)

    # Implement and fit classifier
    clf = classifiers.svm()
    clf.fit(X_train, y_train)
    
    return clf


def evaluation(clf, testing_data, best_k_mers, eval_file):

    # Generate matrices
    best_k_length = len(best_k_mers[0])
    X_test, y_test = matrices.generateXYMatrice(testing_data, best_k_mers, best_k_length)

    # Realize prediction
    y_pred = clf.predict(X_test)

    # Print results
    print("\nClassification report of evaluation\n")
    print(classification_report(y_test, y_pred, digits = 3))
    print("Accuracy :", accuracy_score(y_test, y_pred) * 100, "%\n")
    print("Confusion matrix \n", confusion_matrix(y_test, y_pred))

    # Print to file
    print("Writing to " + eval_file)
    f = open(eval_file, "w")
    f.write("Sequence id, Predicted class \n");
    for i, d in enumerate(testing_data): f.write(d[0] + ", " + y_pred[i] + "\n");
    f.close()

    return y_pred


def prediction(clf, testing_data, best_k_mers, pred_file):

    # Generate matrices
    best_k_length = len(best_k_mers[0])
    X_test = matrices.generateMatrice(testing_data, best_k_mers, best_k_length)

    # Realize prediction
    y_pred = clf.predict(X_test)

    # Print to file
    print("Writing to " + pred_file)
    f = open(pred_file, "w")
    f.write("Sequence id, Predicted class \n");
    for i, d in enumerate(testing_data): f.write(d[0] + ", " + y_pred[i] + "\n");
    f.close()

    return y_pred
