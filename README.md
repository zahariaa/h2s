# hypersphere2sphere (h2s)

hypersphere2sphere(h2s) is a visualization method for conveying the relationships between high-dimensional, labeled data distributions. h2s summarizes each labeled distribution as a sphere (in 3D) or circle (in 2D).

## Basics

The h2s algorithm proceeds in two steps:

#### Step 1: High-dimensional hypersphere estimation
A uniform *N*-ball distribution is fit to the distribution of samples from each category, yielding a center (a *N*-vector) and a radius (a scalar) for each category \cite{Ritter1990,Larsson2008}.

Here `data` is either an [*P* by *N*] matrix or [*P* by *N* by *F*] tensor, where *P* is the number of points (or "samples"), *N* is the dimensionality of the points, and for dynamic data, *F* is the number of frames (or samples in time).

```matlab
data % an [P by N] or [P by N by F] matrix
% create a Hypersphere object of the high-dimensional representation of the data
hi = Hypersphere.estimate(data);
```

The result is a `Hypersphere` object which contains an estimate of the center and radius of the distribution.

```matlab
>> hi
>> hi.categories
```

Running `Hypersphere.estimate` with only one argument assumes that all the data belong to a single label/distribution. You may notice in the `Hypersphere` object, another `Categories` object within it, here in `hi.categories`. This was automatically generated, but you can alternatively specify the labels of each point by first creating your own `Categories` object.

```matlab
cats = Categories(vectors);
```

A minimal `Categories` object simply takes in `vectors`, a [*P* by *C*] logical matrix where *C* is the number of categories, and the (*p*,*c*)th element in the matrix is `true` if the *p*th point belongs to the *c*th category. This matrix is stored in `Categories.vectors`. Labels and colors can be optionally specified in the same object, in `Categories.labels` and `Categories.colors`, respectively.

```matlab
cats = Categories(vectors,labels,colors);
hi = Hypersphere.estimate(data,cats); % cats is now embedded in hi
```

Each of the hyperspheres estimated are in this single `Hypersphere` object. Things are more interesting now, however, with multiple labels. The `Hypersphere` object can lazily compute the statistics of interest that are used to optimize the visualization's appearance in step 2.

```matlab
>> hi
>> hi.overlaps
```

For reference, you can also plot the original points with the colors specified in the `Categories` object. You can also generate new samples from the estimated distributions.

```matlab
% plot the first 2 or 3 dimensions of the original points
hi.plotSamples(points);
% generate arbitrary samples from the hypersphere distributions
hi.plotSamples(100); % 100 samples per hypersphere
```

#### Step 2: Optimization of low-dimensional rendering
If the number of categories *C* is no greater than *n+1* (e.g., for up to 4 categories for a 3d visualization), the hyperspheres' parameters can be perfectly expressed by the visual language of spheres (or circles for a 2d visualization). If *C>n+1*, then a perfect expression of the summary statistics is not in general possible, and the sphere parameters are optimized to best express the hypersphere parameters.

For the second step, H2S initializes the low-dimensional sphere embedding by positioning the sphere centers using MDS with metric stress \cite{Young1938,Torgerson1952,Shepard1962} as the optimization criterion. The sphere embedding configuration is further optimized to minimize the error *E* between the visualized spreads $\tilde{s}$, separation distances $\tilde{d}$, and overlaps, $\tilde{o}$, and the target values measured in the high-dimensional space ($\hat{s}$, $\hat{d}$, and $\hat{o}$, respectively):
<img src="https://latex.codecogs.com/svg.latex?\large&space;E = \sum_{i=1}^{C-1} \sum_{j=i+1}^C \left( \tilde{d}_{ij} - \hat{d}_{ij} \right)^2 + \sum_{i=1}^{C-1} \sum_{j=i+1}^C \left( \tilde{o}_{ij} - \hat{o}_{ij} \right)^2 + \sum_{i=1}^C \left( \tilde{s}_i - \hat{s}_i \right)^2" title="\largeE=\sum_{i=1}^{C-1}\sum_{j=i+1}^C\left(\tilde{d}_{ij}-\hat{d}_{ij}\right)^2+\sum_{i=1}^{C-1}\sum_{j=i+1}^C\left(\tilde{o}_{ij}-\hat{o}_{ij}\right)^2+\sum_{i=1}^C\left(\tilde{s}_i-\hat{s}_i\right)^2" />

This step yields a `SetOfHyps` object (it can be run on a `Hypersphere` or `SetOfHyps` object), which is a subclass of `Hypersphere` that inherits its methods and has additional methods and properties related to h2s optimization, display, and statistical inference. Some properties, such as centers and radii, are protected in `SetOfHyps`, and several other properties (such as overlaps) are eagerly computed and stored in the object.


```matlab
% Run the h2s algorithm to create a low-dimensional representation (converts to a SetOfHyps object)
lo = hi.h2s;       % 3D by default, but...
lo2d = hi.h2s(2);  % ...you can specify a 2D representation

% Render the visualization
lo.show;

% If you so desire, you can do it all in one line:
Hypersphere.estimate(data,cats).show; % (show runs h2s automatically if there are too many dimensions)
```

The H2S visualization may display a black line in the lower right corner. This is an error bar: it represents the maximum error of any of the spread, separation, and overlap terms in the equation above. If H2S arranges the visualization in a way that fails to capture the sign of an overlap/margin correctly, it will place a black line above that overlap to indicate that it is truly a margin, or vice-versa. The `msflips` property identifies which of these overlap/margins have flipped.


## Getting fancy: user options for estimation and optimization

#### Estimation options

The default estimator for the high-dimensional hypersphere parameters is an empirically-derived one which automatically selects the more appropriate of two Minimum Variance Unbiased Estimators (MVUE), one which makes a hyperspherical assumption, and another which makes a Gaussian assumption. If so desired, an alternative estimator can be used, such as one that maximizes the likelihood of a minimally enclosing hypersphere (`'jointml'`), an Markov Chain Monte Carlo (`'mcmc'`) search for the parameters that maximize the posterior probability of the data under the assumption of a uniform *N*-ball distribution, or the MVUEs which make specific distributional assumptions (`'gaussian'`, `'uniformball'`, `'uniformcube'`).

```matlab
hi = Hypersphere.estimate(data,cats,'mcmc');
```

#### H2S optimization options

By default, H2S initializes the low-dimensional sphere embedding by positioning the sphere centers using MDS with metric stress. One can explicitly call for this with the `'mdsinit'` option or, alternatively, use random initializations with `'randinit'`.

```matlab
% runs 10 optimizations with random initializations,
% selects the one with the lowest error.
% the preceding 3 indicates the visualization dimensionality
lo = hi.h2s('randinit',[3 10]);
```

To create an H2S visualization that copies the estimated radii rather than optimizing them jointly with the other summary statistics, one can use regularization terms on a generalized form of the H2S objective:
<img src="https://latex.codecogs.com/svg.latex?\large&space;E = \alpha \sum_{i=1}^{C-1} \sum_{j=i+1}^C \left( \tilde{d}_{ij} - \hat{d}_{ij} \right)^2 + \beta \sum_{i=1}^{C-1} \sum_{j=i+1}^C \left( \tilde{o}_{ij} - \hat{o}_{ij} \right)^2 + \gamma \sum_{i=1}^C \left( \tilde{s}_i - \hat{s}_i \right)^2" title="\largeE=\alpha\sum_{i=1}^{C-1}\sum_{j=i+1}^C\left(\tilde{d}_{ij}-\hat{d}_{ij}\right)^2+\sum_{i=1}^{C-1}\beta\sum_{j=i+1}^C\left(\tilde{o}_{ij}-\hat{o}_{ij}\right)^2+\gamma\sum_{i=1}^C\left(\tilde{s}_i-\hat{s}_i\right)^2" />

```matlab
% disables radius optimization by setting {alpha=1, beta=1, gamma=0}:
lo = hi.h2s({1 1 0});
```

Note that the standard `'mdsinit'` initialization is the global optimum of this reduced objective function.


## Other stuff: statistical inference

The H2S framework can also report statistical significance for each individual summary statistic, and significant differences between them. These are permutation- and bootstrap-based tests that are applied to the high-dimensional data. The output is a `SetOfHyps` object.

```matlab
hi = hi.significance(points); % points is a required input

figure;
ax(1) = subplot(1,3,1);
hi.showSig;           % plots diagram of significance values for each statistic
hi.showSig('legend'); % same as above, but includes a legend indicating thresholded significance values

ax(2) = subplot(1,3,2);
hi.showSig('diff');   % plots diagram of significance of differences among the statistics

% alternatively, this plots both in one line:
hi.showSig(ax);

% you can use the same visualization format to plot the values of the statistics of interest themselves:
ax(3) = subplot(1,3,3);
hi.showValues;

% or plot both the high and low-dimensional values together for comparison
hi.showValues(lo);
```

## Other stuff: dynamics

If you have dynamic data in the form of a [*P* by *N* by *F*] tensor, evaluating H2S on those data, making sure each *f*th frame is optimized in the same visualization space as all the other frames, the code is very much the same:

```matlab
hi = Hypersphere.estimate(data,cats);
lo = hi.h2s;
```

To visualize these frames as frames of a dynamic H2S movie, you can use the `movie` method. The `plotDynamics` method is an alternative visualization that generates static plots of each statistic of interest as they changes from frame to frame.

```matlab
figure;
lo.movie;

figure;
lo.plotDynamics;
```

If you have [*P* by *N* by *F*] tensor data which has the same category structure (i.e., you can apply the same `Categories` object to each frame), but you want each frame to be optimized independently, you can use the `'independent'` option:

```matlab
hi = Hypersphere.estimate(data,cats);
lo = hi.h2s('independent');
```

This can be useful if you've manually bootstrapped or permuted the data.