# Final Project: Mesh Simplification and Remeshing

For more information on this project please visit our personale website: [Google Doc](https://docs.google.com/document/d/1qjw-ya-W3nCAa2ioy9zzPrihgbScVYsUZLmRqXt_J2I/edit).

## Updates

Format: <Date, approximate time - Name>

#### 4/17 5 PM - Kyle
Completed Task 1 and almost done with Task 2. Need to fix bug where vertices with >6 edges get messed up.

Wrote and integrated function definitions for `collapseEdge()`, `downsample()` and `isomesh()`, corresponding to tasks 2, 3 and 4. They can be toggled in the GUI with the hot keys `c`,  `d` and `o`, respectively (`i` was already taken).

#### 4/18 12 PM - Kyle
Fixed bug where collapsing edges connecting to vertices with >6 incident edges would cause holes in the mesh. Edge collapse now appears to be completely functional.
New issue arose where upsampling method appears to not work anymore. Once an edge has been collapsed and you try to upsample it gets stuck in an infinite loop. I imagine this has something to do with the isNew attributes of mesh elements and/or having an irregular mesh structure that causes a problem in the loop subdivision.

#### 4/24 7 PM - Tori
Implemented isotropic remeshing based of of this https://github.com/stanford-cs248/triangle-mesh-editor/wiki/Isotropic-Remeshing explanation. Includes some helper functions, including an implementation for Vertex::computeCentroid() for the Vertex class.
