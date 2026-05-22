/**
 * Solution to https://open.kattis.com/problems/cordonbleu using Hungarian
 * algorithm.
 */

#include <cassert>
#include <iostream>

template <typename T, typename U>
using Pair = std::pair<T, U>;
template <typename T>
using Vector = std::vector<T>;

template <typename T>
using NumericLimits = std::numeric_limits<T>;

/**
 * @brief Checks if b < a
 *
 * Sets a = min(a, b)
 * @param a The first parameter to check
 * @param b The second parameter to check
 * @tparam The type to perform the check on
 * @return true if b < a
 */
template <typename T> 
constexpr bool ckmin(T& a, const T& b) { 
    return b < a ? a = b, true : false; 
}

/**
 * @brief Performs the Hungarian algorithm.
 *
 * Given J jobs and W workers (J <= W), computes the minimum cost to assign each
 * prefix of jobs to distinct workers.
 *
 * @tparam T a type large enough to represent integers on the order of J *
 * max(|C|)
 * @param C a matrix of dimensions JxW such that C[j][w] = cost to assign j-th
 * job to w-th worker (possibly negative)
 *
 * @return a vector of length J, with the j-th entry equaling the minimum cost
 * to assign the first (j+1) jobs to distinct workers
 */
template <typename T> 
Vector<int> hungarian(const Vector<Vector<T>>& C) {
    const int J = static_cast<int>(C.size());
    const int W = static_cast<int>(C[0].size());
    assert(J <= W);
    // job[w] = job assigned to w-th worker, or -1 if no job assigned
    // note: a W-th worker was added for convenience
    Vector<int> job(W + 1, -1);
    Vector<T> ys(J); 
    Vector<T> yt(W + 1);  // potentials
    // -yt[W] will equal the sum of all deltas
    Vector<T> answers;
    const T inf = NumericLimits<T>::max();
    for (int jCur = 0; jCur < J; ++jCur) {  // assign jCur-th job
        int wCur = W;
        job[wCur] = jCur;
        // min reduced cost over edges from Z to worker w
        Vector<T> minTo(W + 1, inf);
        Vector<int> prev(W + 1, -1);  // previous worker on alternating path
        Vector<bool> inZ(W + 1);     // whether worker is in Z
        while (job[wCur] != -1) {    // runs at most jCur + 1 times
            inZ[wCur] = true;
            const int j = job[wCur];
            T delta = inf;
            int wNext;
            for (int w = 0; w < W; ++w) {
                if (!inZ[w]) {
                    if (ckmin(minTo[w], C[j][w] - ys[j] - yt[w]))
                        prev[w] = wCur;
                    if (ckmin(delta, minTo[w])) 
                        wNext = w;
                }
            }
            // delta will always be nonnegative,
            // except possibly during the first time this loop runs
            // if any entries of C[jCur] are negative
            for (int w = 0; w <= W; ++w) {
                if (inZ[w]) {
                    ys[job[w]] += delta;
                    yt[w] -= delta;
                } else {
                    minTo[w] -= delta;
                }
            }
            wCur = wNext;
        }
        // update assignments along alternating path
        for (int w; wCur != W; wCur = w) 
            job[wCur] = job[w = prev[wCur]];
        answers.push_back(-yt[W]);
    }
    return job;
}

/**
 * @brief Performs a sanity check for the Hungarian algorithm.
 *
 * Sanity check: https://en.wikipedia.org/wiki/Hungarian_algorithm#Example
 * First job (5):
 *   clean bathroom: Bob -> 5
 * First + second jobs (9):
 *   clean bathroom: Bob -> 5
 *   sweep floors: Alice -> 4
 * First + second + third jobs (15):
 *   clean bathroom: Alice -> 8
 *   sweep floors: Carol -> 4
 *   wash windows: Bob -> 3
 */
void sanityCheckHungarian() {
    Vector<Vector<int>> costs{{8, 5, 9}, {4, 2, 4}, {7, 3, 8}};
    assert((hungarian(costs) == Vector<int>{5, 9, 15}));
    std::cout << "Sanity check passed." << std::endl;
}
