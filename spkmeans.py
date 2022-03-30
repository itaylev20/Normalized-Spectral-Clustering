import myspkmeanssp
import argparse
import numpy as np


def distance(point1,point2):
    return np.linalg.norm(point1-point2)**2

def select_random_centroid(minimum_distances,n):
    sum_distances=minimum_distances.sum()
    weights= minimum_distances/sum_distances
    return int(np.random.choice(n,p=weights))

def selelct_initial_centroids(points,k,n):

    points_minimum_distances = np.full(n, np.inf)
    current_centroid_index=int(np.random.choice(n))
    centroids_indices = [current_centroid_index]

    for z in range(2,k+1):

        for index, point in enumerate(points):

            distance_from_last_centroid=distance(point,points[current_centroid_index])
            if points_minimum_distances[index] > distance_from_last_centroid:
                points_minimum_distances[index] = distance_from_last_centroid

        current_centroid_index=select_random_centroid(points_minimum_distances,n)
        centroids_indices.append(current_centroid_index)


    return centroids_indices


def printVectors(vects):
    for v in vects:
        print(','.join([str(np.round(x,decimals=4)) for x in v]))

goals={'spk':None,'wam':None,'ddg':None,'lnorm':None,'jacobi':None}
max_iter=300

if __name__ == '__main__':
    np.random.seed(0)
    parser = argparse.ArgumentParser(description='Spkmeans')
    parser.add_argument('k', type=int, help="Clusters number")
    parser.add_argument('goal', type=str, choices=goals.keys() , help=f"Can get the following values {goals.keys()}")
    parser.add_argument('file_name', type=str, help="The path to the file that will contain N observations, the file extension is .txt or .csv.")
    args = parser.parse_args()

    k = args.k
    file_name = args.file_name
    goal = args.goal

    points=np.loadtxt(file_name, delimiter=",",ndmin=2)
    n=len(points)
    d=len(points[0])

    assert (0 <= k < n), "Invalid input."

    # argparse take care of invalid choice of goals - therefore we can be sure the string goal is valid
    if goal=='spk':
        T,k =myspkmeanssp.getT(k,n,d,points.tolist())
        T=np.array(T)
        initialCentroidsIndicesNumpyArray = selelct_initial_centroids(T, k, n)
        initialCentroids = T[initialCentroidsIndicesNumpyArray]

        print(','.join( str(int(x)) for x in initialCentroidsIndicesNumpyArray))
        myspkmeanssp.kmeans(k,n,k,max_iter,initialCentroids.tolist(),T.tolist())

    else:
        myspkmeanssp.sideGoals(n,d,goal,points.tolist())


# todo - both python and C - output last line should not be with \n Discussion https://moodle.tau.ac.il/mod/forum/discuss.php?d=161522
# Todo - does eigenvectors are normalized? did https://moodle.tau.ac.il/mod/forum/discuss.php?d=161539יש
# TODO https://moodle.tau.ac.il/mod/forum/discuss.php?d=162884 ?
# TODO -  jacobi -  open discussion about last iteration - https://moodle.tau.ac.il/mod/forum/discuss.php?d=162884
# TODO - How rami will know which files to compile ? https://moodle.tau.ac.il/mod/forum/discuss.php?d=156940
# TODO - open discussion about %4f and rounding  https://moodle.tau.ac.il/mod/forum/discuss.php?d=162710