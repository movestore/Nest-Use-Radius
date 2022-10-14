# Nest Use by Radius
MoveApps

Github repository: *github.com/movestore/Nest-Use-Radius*

## Description
Calculates the proportion of locations/duration of your track in one (or more) given radius/radii around a nest/central location. The latter have to be added to the data in a previous App. Also TimeLag has to be added to the data before.

## Documentation
This App takes as input a nest/central location for each track and a list of radii (provided by the user). Then, the proportion of locations and duration (time lag see below) of each track in each provided radius around the nest/central location is determined and returned in a csv table. For each radius also mean and standard deviation are provided.

In addition, for each track, a map is returned indcluding the nest location and circles indication the different radii around it.

For durations, one has to include the App "Time Lag Between Locations" in the workflow, best before filtering by speed (as otherwise the duration become unrealistically too long).

Note that the nest/central location is taken as the coordinates of the respective object.

### Input data
moveStack in Movebank format

### Output data
moveStack in Movebank format

### Artefacts
`Radius_NestUse.csv`: table of proptions of locations/durations for each individual and radius. Per radius also averages (mean) and standard deviation are given.

`Tracks_withRadii_onMap.pdf`: For each track an openstreet map with the track plotted together with the nest/central location and radii indicated by circles around it.

### Parameters 
`radii`: One or more radius/radii (in metre) that you want proportional use calculated for. For multiple values please separate by comma. Default 500.

`selVar`: Name of the object containing the nest/central location coordinates. E.g. `nesting` when using csv output from Nest Detection App.

`trackVar`: Name of the track ID variable in the object containing the nest/central location coordinates. Take care that this parameter also exists in the track attributes of the input data set.


### Null or error handling:
**Parameter `radii`:** This parameter must contain positive values only. If multiple radii are not separated by comma, an error will occur.

**Parameter `selVar`:** If there is no object with the given name the App will run into an error.

**Parameter `trackVar`:** If there is no variable with the name given here, an error will be returned.

**Data:** The full data set is returned.