## CS184 Final Project 
# Mesh Simplification and Remeshing
Victoria Austin, Kyle Rentschler, Max Yao

# Summary 
In this project we aim to build on the topics covered in assignment 2 and add new functionalities to our mesh editor, specifically downsampling using quadric error simplification and isotropic remeshing. Downsampling reduces the number of polygons in a mesh which results in less storage being used, faster/cheaper rendering and simpler manipulation of objects. Remeshing is the process of refactoring a mesh structure such that the number of triangles remains the same but the size/shape of each triangle is more regular/standardized. These two methods can help us build more accurate and simpler mesh structures for use in rendering applications.

# Problem Description
Downsampling a mesh is the process of removing mesh features (edges, faces, vertices) to result in a simplified mesh structure using fewer triangles. Reducing the number of polygons to represent an object is very useful for when you want to trade detail for speed, such as rendering an object that is far away or small. Removing triangles requires the use of basic mesh operators, such as an edge collapse method. When an edge is collapsed, two vertices (and some of their connecting edges) are joined into one. This creates a problem of how to decide where to place the resulting merged vertex in 3D space. Another issue that downsampling introduces is deciding which edges to collapse in a way that preserves the original structure and keeps as many details as possible. In other words, we want to remove triangles that are the most redundant. This is a difficult task because the problem space is massive; at any step in the downsampling process you can remove any edge and place the resulting vertex anywhere, which all begins again when you have to remove the next edge. Error quadrics solves these problems by allowing us to quantify the “cost” of removing an edge, as well as determine the ideal location to place the vertex resulting from the collapse such that it minimizes distance to the original surface. Though it is only an approximation, it is a very good approximation and, unlike finding the true optimum, we can solve for this efficiently.

Another factor to consider is that downsampling can lead to ill-formed meshes with highly-irregular triangle structures. This can also occur when meshes are initially formed, such as when they are obtained through a 3D scan. This is where remeshing comes into play. Remeshing is the process of refactoring a mesh structure such that the number of triangles remains the same but the size/shape of each triangle is more regular/standardized. Note that “improved” is somewhat arbitrary here, and exploring what is one of the questions we hope to address in this project. The main tools at our disposal to address this issue are edge manipulation (flips, deletes, splits) and vertex manipulation (moving the position of the vertex). This task is challenging because, like downsampling, the output space for all these changes is infinite and there is no strong definition of what improvement is. Good heuristics are needed to use these tools well in combination, which if done can result in a much more regularized mesh with little deviation from the original shape. Our approach to solving this problem will be to research, test and implement an effective algorithm (including heuristics) for how to apply the remeshing fundamentals. Possibilities for our solution include enforcing a maximum/minimum edge length (and splitting or collapsing edges that exceed this threshold), rules for repositioning vertices and flipping edges. 

# Goals and Deliverables:
What We Plan to Deliver:
Because we are working with mesh structures and the nature of our improvements are subjective with no quantitative way to measure improvement, our base deliverables will be in the form of a working code base implementing these techniques and sample images showing our results. Our goal deliverables include developing a way to rate mesh structure changes numerically to include in our analysis. 

Task 1: Code reformatting
Our first task will be to modify the initial . The following code deliverables will need to be completed:
A function definition and corresponding hotkey to collapse an edge (use edge split as a model for how this can be used with a selected edge in the GUI with a key press and also be called by other parts of the code to do programatically when we implement downsampling)
A function definition and corresponding hotkey to run an iteration of downsampling via quadric error metrics.  Possible implementation additions could include a parameter to run n iterations of this algorithm each step. 
A function definition and corresponding hotkey to run the remeshing algorithm. Possible implementation additions could include some kind of parameter to adjust heuristics (e.g. specifying tolerance for edge length variance, the number of resulting triangles, etc).

Task 2: Edge collapse
Implement edge collapse in the corresponding function definition provided from task 1. Make sure to include contingencies for edge cases (ignore requests to collapse boundary edges) and test sufficiently to make sure no bugs are present because this function will be heavily relied on in the later parts. For now, just assign the position of the new vertex to the midpoint of the edge being collapsed.
The following deliverables will need to be completed:
Filling in the function definition from task 1 with logic to implement edge collapse
A series of at least two pairs of before/after images of a mesh showing an edge being collapsed

Task 3: Downsampling via Quadric Error Metrics
Implement downsampling via quadric error metrics. The first time this is called it will need to create the priority queue of the order of edge priority (where priority is determined by quadric error). Each time it is called it will pop off the top item from the queue and collapse the edge, and then update the queue. Possible extensions include the option to do this in batches (i.e. multiple edges at a time) and also including calls to the remeshing function periodically throughout the downsampling process.
Deliverables:
Filling in the function definition from task 1 with logic to implement quadric error metrics and collapse one edge.
A series of at least three pairs of before/after images showing the results from sustained downsampling (i.e. many triangles reduced) of different meshes.

Task 4: Isotropic Remeshing
Implement the isotropic remeshing function. Make sure to include edge splitting, edge collapsing, edge flipping and vertex repositioning when appropriate. Possible extensions include the possibility to adjust heuristics according to user-given parameters.
Deliverables:
Filling in the function definition from task 1 with logic to implement isotropic remeshing and fully remesh and entire mesh structure in one call.
A series of at least three pairs of before/after images showing the results from remeshing different meshes.

# Questions for analysis:
What is the relative tradeoff of number of triangles to level of detail? E.g. is there a significant dropoff in quality after a certain percentage reduction in complexity?
When is it advantageous to standardise the shape and size of our triangles with remeshing and when would it be better to leave it be? 
Is it helpful to the final result if we do intermediate remeshing while doing quadric error downsampling
How do different features or shapes cope with downsampling? 

# What We Hope to Deliver:
If everything runs smoothly, we anticipate being able to implement one of the three tasks described below. We hope to implement one of these tasks in order to provide further comparisons between different downsampling algorithms and to be able to discuss the tradeoffs between the efficiency, standardization, and geometry/feature preservation of downsampling a mesh structure. 

Task 1: Vertex Decimation: 
Vertex decimation is another downsampling algorithm but is slightly more complex than quadric error simplification, since it has more edge cases to consider. 
Implement vertex removal, which is similar to edge collapse but focuses on removing vertices rather than edges. To do this, we would have to create a vertex classification scheme that characterizes which vertices are the most advantageous to remove and place them in a priority que for removal. 
Implement re-triangulation, which fills holes formed from removing vertices.
This algorithm has the possibility of changing the topology of the mesh and making it non-manifold. To avoid this, we must provide additional checks that ensure that we do not create duplicate triangles and edges.  
A comparison between our two downsampling algorithms, including run time analysis and overall differences in the quality and structure of our downsampled meshes through a series of before/after images.

Task 2: Geometric Error Quantifier: 
We would like to implement a quantitative metric that describes our algorithms' effectiveness at producing a downsampled mesh that is as close to the original geometry of the mesh as possible. One possible method for this is to compute the Hausdorff distance between the original mesh and the downsampled mesh. The Hausdorff distance measures the difference between two different representations of the same 3D object by finding the greatest distance from a point in one mesh to the closest point in the other mesh. 
Implement finding the Hausdorff distance between two meshes. 
Compare how using different methods for downsampling affect the Hausdorff distance between two meshes, and discuss the possible trade offs between efficiency and preserving the geometry across our various downsampling methods 

Task 3: Fast Mesh Simplification with Feature-Preserving Efficiency :
We would like to explore the tradeoffs between standardizing the shape and size of triangles in a mesh using Isotropic remeshing and using another Feature-Preserving algorithm that omits excessive details and preserves salient features within a mesh during downsampling. 
Implement the calculation of the Energy-Operator 
Implement triangle decimation and new vertex calculation using our Energy-Operator heuristic 
Implement topological reconstruction to spatially connect new and old vertices to preserve the topology of regular manifold meshes 
 A comparison between our two downsampling algorithms, including run time analysis and overall differences in the quality and structure of our downsampled meshes through a series of before/after images. 

# Schedule:
Proposal Due: Thursday April 9th 11:59pm.
Proposal Revisions Sent, Project Start: Tuesday April 14th.
Graded Milestone Due: Tuesday April 28th 11:59pm.
Final Presentations: Thursday May 7th.
Final Deliverables Due: Tuesday May 12th 11:59pm
Week 1: (4/14) - (4/21)
Implement any revisions suggested by course staff
Complete task 1
Complete task 2
Begin researching and write pseudocode for tasks 3 and 4
Week 2: (4/21 - 4/28)
Do bulk of work for tasks 3 and 4
Begin collecting sample images for report and writing about results
Debug time for tasks 3 and 4
Seek out course staff support (if necessary)
Complete and submit milestone deliverables

Week 3: (4/28 - 5/5)
Complete tasks 3 and 4
Generate sample images for downsampling and remeshing for report
Practice final presentation
Write up analysis of tasks and results
Begin work on reach goals (time permitting)

Week 4: (5/5 - 5/12)
Prepare all deliverables for final report
Continue working on reach goals (time permitting)
Write about reach goals (if achieved)

# Resources: 
The primary inspiration for this project is from this assignment from Carnegie Mellon University. For software resources, we will be starting with a completed CS184 Assignment 2 codebase from one of our projects and utilizing all the same libraries and frameworks in the code base along with the sample mesh files and GUI. The computing platform and hardware we will use are personal machines. In addition, we will be using the following online resources (more to potentially follow).
This paper on quadric error metrics by Garland and Heckbert, the original creators of this method
This lecture from Stanford on isotropic remeshing 
This lecture from USC on isotropic remeshing 
This lecture from Princeton to help with quadric error metric implementation
Possibly use this lecture from Michigan Tech on vertex decimation if we decide to pivot or include it as a reach goal 
Possibly us this https://tel.archives-ouvertes.fr/tel-00838783/file/these.pdf paper on how to calculate the Hausdorff distance between two meshes 
Possibly use this http://downloads.hindawi.com/journals/sp/2019/4926190.pdf paper on Fast Mesh Simplification with Feature-Preserving Efficiency 

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
