# Background Concentration
One of the fundamental assumptions of the GIM method is that we can describe the observed methane concentration as the sum of three factors: emissions propagated through the atmosphere (from sources both in and out of the scene), an underlying latent background concentration of methane, and instrument noise. In order to run our inversion methods, we need to somehow separate these factors, which means coming up with an estimate for the background, ideally with some measure of uncertainty to constrain the final estimate. Since spatial misallocation of the background concentration can propogate to retrieved emissions, it is important to fit the proper magnitude and spatial distribution of the background. The current leading theory relies on our understanding of these factors:

1. **Emissions**: Emissions are strictly positive rates of flux that blow thorugh our scene. We make some heavy assumptions about the emissions for our inversion - that their motion is perfectly characterized by the Jacobian and that they are constant in time. For the purposes of Background Estimation, we make less ambitious assumptions - we assume that emissions are strictly positive, and that there are one or more measurable areas of every scene that have no emissions. The latter assumption might need some justification - it is both by design of the MSat targets and by experience that we claim this to be true. 

2. **Background**: The concept of a background is an artificial construct designed to capture any gradients not attributable to emissions in the Total (Jacobian) Domain. This is like asking - if there were no emissions, what would methane concentrations look like? Our theory is that the background concentration is equal to the L2 Prior plus a scene-specific constant, modulated by the averaging kernel (the averaging kernel is to account for instrument sensititvity). The L2 Prior is a model that accounts for various factors - including topography, latitude, date, tropopause height - to give a met-informed but data-blind estimate. L2 uses it as a prior for a Bayesian update process to fit to the observed spectra. While the averaging kernel varies spatially and vertically, we make the simplification that the prior offset from true values happens equally in all atmosphere levels. Importantly for us, this means that we need to only fit a single value in our Background Estimation process - the additive offset, *c*.

```math
\mathbf{XCH4}_{bg} = \mathbf{XCH4}_{prior}  + c \mathbf{A} + \epsilon
```

where bold values denote matrices (image arrays), ***A*** is the column-averaged averaging kernel, and the true XCH4 in background regions is *c* ppb away from the L2 Prior for every element *i* in the image : $`c = \mathbf{XCH4}_{true, i} - \mathbf{XCH4}_{prior, i}\ ;    \ \ \forall i \in \mathbf{XCH4}_{prior}`$

3. **Noise**: The last factor, ε, is noise, whether from the instrument, retrieval forward model, imprecision of our math, incorrect retrieval parameters like surface pressure and degree of specular reflection, aerosols, cloud 3D effects, atmospheric mixing, or some other factor. We believe the instrument noise to be normal, centered at 0, and that its scale can be broadly predicted by the albedo of the scene. This assumption likely doesn't hold for all components of the noise. For an example noise model for TROPOMI, see Hu et al. (2016).

<img src="../doc_resources/Instrument-noise-signal.png" width="600">

Instrument noise as a function of radiance [which scales with albedo * cos(solar zenith angle)]. This figure was made from 268 L3s, filtered by O2DP and total pixel count, after data domain selection. Error bars show the standard deviation of each bin.

## Fit from below (FFB) Methods

Taken together, this leads us to a hypothesis: there is some part of the scene, where emissions do not have a sizable effect, where observed concentrations should look like a normal centered on the L2 Prior plus some constant `c` times the averaging kernel. If we could find those areas, we could take a mean and subtract the Prior, divide by the averaging kernel, and get `c`. 

The *Reflected Normal Noise* (RNN) approach is an attempt to do so. We note that if we knew the correct background value, everything that we observe below that value would be noise. Assuming noise is normal with a predictable scale, we try to find a background value such that all values below it (the residuals) best fit this model. Other methods fit the expected standard deviation, either derived from a noise model or empirically.

The background concentration is fit by:
1. Calculating the retrieval error for each observation
2. Selecting a range of possible values for the offset, *c*
3. For each possible value of *c*, computing the candidate background from the equation above and its reflected normal noise
4. Selecting the value of *c* that best passes a noise test

For convenience, the methods in `negative_reflection.py` take a given scene and return a distribution of data that represents the full normal, by taking all values **below** the proposed background and then reflecting across the `x=0` vertical. 

### Normality Tests

One simple approach is to test this distribution with established statistical metrics of normality, and to pick the background that leaves the most normal-esque negative reflection. We use Scipy's built in `normaltest` to do this, because it's fast and looks pretty identical in performance to other normality tests. 

### Z-sigma normalization tests

Another approach is to compare the spread of the negative residuals (i.e., the observations minus a proposed background) to the spread we expect from the observation precision. Here, spread means standard deviation. Ideally, negative residuals are only due to noise in the data. Deriving the expected spread for the data from first principles (instrument noise, ...) has proven prohibitively hard. Therefore, we estimate the instrument noise from local neighborhoods in the data (given a constant, true value, the spread of the observation is the noise). We pick the background proposal that produces a negative residual distribution that best matches our expectations, i.e., the negative deviations divided by the local standard deviation have a spread of one.
Given some L3 xCH4 observations and the mean averaging kernel at those points, estimate the background using the Z-Sigma method performs these steps:

1. Calculate the local std and mean for the whole scene at the L3-resolution using a local neighborhood, e.g. a 2 km x 2 km square.
2. For each background guess, take the negative residuals, divide them by the local std, and reflect them around 0 (i.e. RNN normalized to local std).
3. Calculate the std of the RNN. The best background is the one with the std closest to 1.

The background may, depending on the data, have an obvious bias. Connected patches of negative residuals are a quality criterion for the background estimate. If the largest cluster is very small, the background is probably underestimated. If the largest cluster is very large, the background is probably overestimated. We can try to salvage the background estimate by iteratively optimizing the background guess:

1. Large clusters: Probably due to enhanced regions contributing to negative residuals. -> Iteratively remove data in regions with high local means. An area with an enhanced local mean is less likely to be a background region, but can have negative residuals due to noise. Since we remove data from high-mean regions, this should be fine.
2. Small clusters: Probably due to low-lying outliers in the scene. -> Iteratively remove the smallest values from the enhancements. This is very sketchy, and I have much less confidence in this approach than in the large cluster approach since we operate explicitly on the lower end of the data distribution.

### Regional Tests

Some tests, notably the `normality` test, have a regional variant. This takes spatially connected chunks of the L3 data and runs the method on them, and then takes a 20th percentile (arbitrary number) among the areas.

# Code Structure

The only entrypoint for production code is in `background/api.py` and is through the `estimate_background_for_l3` function. This function will in turn parse the options passed to it to choose a specific method to run, and execute that method directly. 

## Adding a New Method

If you are adding a new method, the intended workflow is to create a new enum for it in `background/data_types.py` in `BackgroundCalculationMethod`. Then, if your method requires new options, modify `BackgroundCalculationOptions` - note that new values in this dataclass should have default values and should likely be `Optional` as long as it is only appropriate for some methods. 

You will likely want to create a new file in this directory `background/cool_method.py` and at least one function that follows the signature of other method functions - it takes the L3 regrid data offset from the prior, the averaging kernel, and an options struct; it outputs `BackgroundCalculationResults` and `BackgroundCalculationAnalytics`. For convenience, if your method is a statistical test of some kind, applied either on the whole scene or convolutionally, you may benefit from calling the `whole_scene_test` and `convolutional_test` functions from `background/utility.py`, which handle many common operations. You will still have to aggregate your convolutional results, however. 

### References

Hu, H., Hasekamp, O., Butz, A., Galli, A., Landgraf, J., Aan De Brugh, J., Borsdorff, T., Scheepmaker, R., & Aben, I. (2016). The operational methane retrieval algorithm for TROPOMI. Atmospheric Measurement Techniques, 9(11), 5423–5440. https://doi.org/10.5194/amt-9-5423-2016
