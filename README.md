# Morphometrics wrapper package for consistent analyses of 3D Shape. - WIP
Uses Morpho's and geomorph's core functions. Contains several basic functions that can be used as is and larger wrapper functions that connect these basic functions to an upcoming UI interface.
This package is a work-in-progress and has only been tested on Ubuntu.

**Basic functions**
- convert() Detects format of shapes and converts between 1d <-> 2d <-> 3d
- residualise() Given select covariates, uses evaluate_model() to find best fitting linear model and then corrects shapes according to that model (residuals is added back to the mean)
  - evaluate_model() Given a linear model, detect which covariates are significant to find best fitting linear model
- scrutinise() Test for differences in shape (or size) by independent variable(s) with a Procrustes ANOVA (geomorph::procD.lm). If there is one independent variable and the model is significant: post-hoc pairwise testing (if >2 levels) and CVA morphs (categorical), or morphs along the variable at +/- 3SD (continuous).
- trim() Reduce dimensionality of shapes by parallel analysis (PCAtools::parallelPCA, default), scree plot (interactive), or a cumulative variance threshold
- compare_axes() Measure alignment between mean shapes. - WORK IN PROGRESS
- caliper() Measure linear distances in mouse crania with MusPO::MeasureAll(). The group is the EOS independent variable. Requires derive_meta() only, and can be run at any later stage. Only functional using BHLab's landmark set.

**Visualisation functions**
*These are ggplot2 bassed figures*
- get_palette() Vector of colours associated with default plots from this package
- get_theme() ggplot2 theme - mask this theme if you have your preferred theme
- plotbars() bar plot - gpglot2 object
- plotdensity() density plot - ggplot2 object
- plotscatter() scatter plot - ggplot2 object
- biplot() scatter plot + x and y density plots - ggarrange object
*Morphing functions*
- shapeshifter() Given landmarks and a mesh (e.g. atlas), deforms mesh to landmarks (can take meshes with or without the same # of vertices as landmarks)
- heatmapify() Create heatmap morphs using Blue-Red or Green-Magenta color palettes

**Specific functions for BHLab**
- query_data() Query data from BHLab servers. Set server location using set_servers(). For mouse data, data is read using import_mm_data(). For HS data, it is read from prepared "Basic" files. This function will be overhauled when storages are restructured.
  - import_mm_data() Read mouse data (.tag files) directly from servers - this is an edit of Marta Vidal Garcia's tag2array.
  - set_servers() Set server directory on local computer.

**Functions for UI interface**
*These functions operate based on an "EOS" class of object*
- initialise_eos() Format landmarks and metadata into an EOS object; pass parameters for analysis.
- derive_meta() Detect metadata columns (e.g. Sex, age, etc.); Removes constant variables from meta data. Indicate dependent and independent variables here.
- process_shapes() Symmetrise (Morpho::, align (Morpho::procSym), and reduce dimensions of shape (reduce_shapes). Shape effects of select covariates as indicated in derive_meta() are removed using residualise().
- find_diffs() Determine differences in size and shape by looping through scrutinise(). Uses procD.lm for AOV. If there is one independent variable: conducts post-hoc testing (if multiple levels). If shape effects are significant, then how groups differ will also be examined (e.g. CVA  or linear model).
- subsetting method: [subjects,landmarks]
