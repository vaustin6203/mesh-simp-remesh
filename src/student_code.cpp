#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  {
    // TODO Part 1.
    std::vector<Vector2D> returnPoints;
    for (int i = 0; i < points.size() - 1; i++) {
        Vector2D intermed = (1-t)*points[i] + t*points[i+1];
        returnPoints.push_back(intermed);
    }
//    return std::vector<Vector2D>();
    return returnPoints;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> returnPoints;
      for (int i = 0; i < points.size() - 1; i++) {
          Vector3D intermed = (1-t)*points[i] + t*points[i+1];
          returnPoints.push_back(intermed);
      }
      return returnPoints;
//    return std::vector<Vector3D>();
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> returnPoints = points;
      for (int i = 0; i < points.size() - 1; i++){
          returnPoints = evaluateStep(returnPoints, t);
      }
    return returnPoints[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      std::vector<Vector3D> collapsedPoints;
      for (int i = 0; i < controlPoints.size(); i++) {
          Vector3D intermed = evaluate1D(controlPoints[i], u);
          collapsedPoints.push_back(intermed);
      }
      return evaluate1D(collapsedPoints, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      
      HalfedgeCIter start = this->halfedge();
      
      HalfedgeCIter current = this->halfedge();
      
      std::vector<Vector3D> normalVecs;
      
      do {
          Vector3D v1 = current->vertex()->position;
          Vector3D v2 = current->next()->vertex()->position;
          Vector3D v3 = current->next()->next()->vertex()->position;
          Vector3D crossProd = cross(v2-v1, v3-v1);
          normalVecs.push_back(crossProd);
          current = current->twin()->next();
      } while (current != start);
      int numVecs = normalVecs.size();
      Vector3D avgSum = normalVecs[0] / numVecs;
      for (int i = 1; i < numVecs; i++) {
          avgSum += normalVecs[0] / numVecs;
      }
      return avgSum.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
      if (e0->halfedge()->isBoundary()) {
          return e0;
      }
      
      //halfedge elements
      HalfedgeIter h0 = e0->halfedge();
      HalfedgeIter h1 = h0->next();
      HalfedgeIter h2 = h1->next();
      HalfedgeIter h3 = h0->twin();
      HalfedgeIter h4 = h3->next();
      HalfedgeIter h5 = h4->next();
      HalfedgeIter h6 = h1->twin();
      HalfedgeIter h7 = h2->twin();
      HalfedgeIter h8 = h4->twin();
      HalfedgeIter h9 = h5->twin();
      
      //vertex elements
      VertexIter v0 = h0->vertex();
      VertexIter v1 = h3->vertex();
      VertexIter v2 = h6->vertex();
      VertexIter v3 = h8->vertex();
      
      //edge elements
      EdgeIter e1 = h1->edge();
      EdgeIter e2 = h2->edge();
      EdgeIter e3 = h4->edge();
      EdgeIter e4 = h5->edge();
      
      //face elements
      FaceIter f0 = h0->face();
      FaceIter f1 = h3->face();
      
      //reassign halfedges
      h0->setNeighbors(h1, h3, v3, e0, f0);
      
      h1->setNeighbors(h2, h7, v2, e2, f0);
      
      h2->setNeighbors(h0, h8, v0, e3, f0);
      
      h3->setNeighbors(h4, h0, v2, e0, f1);
      
      h4->setNeighbors(h5, h9, v3, e4, f1);
      
      h5->setNeighbors(h3, h6, v1, e1, f1);
      
      h6->setNeighbors(h6->next(), h5, v2, e1, h6->face());
      
      h7->setNeighbors(h7->next(), h1, v0, e2, h7->face());
      
      h8->setNeighbors(h8->next(), h2, v3, e3, h8->face());
      
      h9->setNeighbors(h9->next(), h4, v1, e4, h9->face());
      
      //reassign vertices
      v0->halfedge() = h2;
      v1->halfedge() = h5;
      v2->halfedge() = h3;
      v3->halfedge() = h0;
      
      //reassign edges
      e0->halfedge() = h0;
      e1->halfedge() = h5;
      e2->halfedge() = h1;
      e3->halfedge() = h2;
      e4->halfedge() = h4;
      
      //reassign faces
      f0->halfedge() = h0;
      f1->halfedge() = h3;
      
//    return EdgeIter();
      return e0;
  }



    VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
    {
        // TODO Part 5.
        // This method should split the given edge and return an iterator to the newly inserted vertex.
        // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

        if (e0->halfedge()->isBoundary()) {
            return e0->halfedge()->vertex();
        }

        //halfedge elements
        HalfedgeIter h0 = e0->halfedge();
        HalfedgeIter h1 = h0->next();
        HalfedgeIter h2 = h1->next();
        HalfedgeIter h3 = h0->twin();
        HalfedgeIter h4 = h3->next();
        HalfedgeIter h5 = h4->next();
        HalfedgeIter h6 = h1->twin();
        HalfedgeIter h7 = h2->twin();
        HalfedgeIter h8 = h4->twin();
        HalfedgeIter h9 = h5->twin();

        //vertex elements
        VertexIter v0 = h0->vertex();
        VertexIter v1 = h3->vertex();
        VertexIter v2 = h6->vertex();
        VertexIter v3 = h8->vertex();

        //edge elements
        EdgeIter e1 = h1->edge();
        EdgeIter e2 = h2->edge();
        EdgeIter e3 = h4->edge();
        EdgeIter e4 = h5->edge();

        //face elements
        FaceIter f0 = h0->face();
        FaceIter f1 = h3->face();

        //declare new halfedges
        HalfedgeIter h10 = newHalfedge();
        HalfedgeIter h11 = newHalfedge();
        HalfedgeIter h12 = newHalfedge();
        HalfedgeIter h13 = newHalfedge();
        HalfedgeIter h14 = newHalfedge();
        HalfedgeIter h15 = newHalfedge();

        //declare new faces
        FaceIter f2 = newFace();
        FaceIter f3 = newFace();

        //declare new vertex
        VertexIter m = newVertex();

        //declare new edges
        EdgeIter e5 = newEdge();
        EdgeIter e6 = newEdge();
        EdgeIter e7 = newEdge();

        //assign vertex position to midpoint
        m->position = 0.5*(v1->position + v0->position);

        //reassign halfedges
        // quadrant 1
        h0->setNeighbors(h1, h3, m, e0, f0);
        h1->setNeighbors(h12, h6, v1, e1, f0);
        h12->setNeighbors(h0, h13, v2, e6, f0);
        h6->setNeighbors(h6->next(), h1, v2, e1, h6->face());

        //quadrant 2
        h2->setNeighbors(h10, h7, v2, e2, f3);
        h10->setNeighbors(h13, h11, v0, e7, f3);
        h13->setNeighbors(h2, h12, m, e6, f3);
        h7->setNeighbors(h7->next(), h2, v0, e2, h7->face());

        //quadrant 3
        h11->setNeighbors(h4, h10, m, e7, f2);
        h4->setNeighbors(h14, h8, v0, e3, f2);
        h14->setNeighbors(h11, h15, v3, e5, f2);
        h8->setNeighbors(h8->next(), h4, v3, e3, h8->face());

        // next, twin, vertex, edge, face
        //quadrant4
        h3->setNeighbors(h15, h0, v1, e0, f1);
        h15->setNeighbors(h5, h14, m, e5, f1);
        h5->setNeighbors(h3, h9, v3, e4, f1);
        h9->setNeighbors(h9->next(), h5, v1, e4, h9->face());

        //reassign vertices
        v0->halfedge() = h4;
        v1->halfedge() = h1;
        v2->halfedge() = h2;
        v3->halfedge() = h5;
        m->halfedge() = h0;

        //reassign edges
        e0->halfedge() = h0;
        e1->halfedge() = h1;
        e2->halfedge() = h2;
        e3->halfedge() = h4;
        e4->halfedge() = h5;
        e5->halfedge() = h15;
        e6->halfedge() = h12;
        e7->halfedge() = h10;

        //reassign faces
        f0->halfedge() = h1;
        f3->halfedge() = h2;
        f2->halfedge() = h4;
        f1->halfedge() = h3;

        //assign isNew values
        m->newPosition = e0->newPosition;
        m->isNew = true;
        e5->isNew = true;
        e6->isNew = true;
        e7->isNew = false;
        e0->isNew = false;


        return m;
    }


//    VertexIter HalfedgeMesh::collapseEdge( EdgeIter e0 )
//  {
//    // This method was added as part of our final project modifications.
//    // It collapses the given edge and return an iterator to the newly inserted vertex.
//
//      if (e0->halfedge()->isBoundary()) {
//          return e0->halfedge()->vertex();
//      }
//
//      //halfedge elements
//      HalfedgeIter h0 = e0->halfedge();
//      HalfedgeIter h1 = h0->next();
//      HalfedgeIter h2 = h1->next();
//      HalfedgeIter h3 = h0->twin();
//      HalfedgeIter h4 = h3->next();
//      HalfedgeIter h5 = h4->next();
//      HalfedgeIter h6 = h1->twin();
//      HalfedgeIter h7 = h2->twin();
//      HalfedgeIter h8 = h4->twin();
//      HalfedgeIter h9 = h5->twin();
//      HalfedgeIter h15 = h6->next();
//      HalfedgeIter h14 = h15->twin()->next();
//      HalfedgeIter h16 = h14->twin()->next();
//
//      //vertex elements
//      VertexIter v1 = h0->vertex();
//      VertexIter v0 = h3->vertex();
//
//      //edge elements
//      EdgeIter e1 = h1->edge();
//      EdgeIter e2 = h2->edge();
//      EdgeIter e3 = h4->edge();
//      EdgeIter e4 = h5->edge();
//
//      //face elements
//      FaceIter f0 = h0->face();
//      FaceIter f1 = h3->face();
//      FaceIter f2 = h7->face();
//      FaceIter f3 = h6->face();
//      FaceIter f4 = h8->face();
//      FaceIter f5 = h9->face();
//
//
//      //declare new vertex
//      VertexIter m = newVertex();
//
//
//      //assign vertex position to midpoint
//      m->position = 0.5*(v1->position + v0->position);
//
//      // vertices
////      h14->vertex() = v1;
////      h15->vertex() = v1;
////      h9->vertex() = v1;
////      h16->vertex() = v1;
//
//      HalfedgeIter curr = h15;
//      do {
//          curr->vertex() = v1;
//          curr = curr->twin()->next();
//      }
//      while(curr != h9);
//      h9->vertex() = v1;
//
//      // twins
//      h6->twin() = h7;
//      h7->twin() = h6;
//      h8->twin() = h9;
//      h9->twin() = h8;
//
//      // edges
//      h6->edge() = e2;
//      h9->edge() = e3;
//
//      // halfedges
//      e2->halfedge() = h6;
//      e3->halfedge() = h9;
//
//      // delete all components that get removed from the edge collapse
//      deleteHalfedge(h0);
//      deleteHalfedge(h1);
//      deleteHalfedge(h2);
//      deleteHalfedge(h3);
//      deleteHalfedge(h4);
//      deleteHalfedge(h5);
//      deleteEdge(e0);
//      deleteEdge(e1);
//      deleteEdge(e4);
//      deleteFace(f0);
//      deleteFace(f1);
//      deleteVertex(v0);
//
//      v1->position = m->position;
//
//
//      return v1;
//  }

    VertexIter HalfedgeMesh::collapseEdge( EdgeIter e0 )
    {
        /**
              v02                v02
               .                 .
              . .                |
             .   .               |
          v00.-----.v01  ->        v00
             .   .               |
              . .                |
               .                 .
               v10                v10
    **/

      /// Phase I: Collect Elements
      // face 0
      HalfedgeIter h00 = e0->halfedge();
      HalfedgeIter h01 = h00->next();
      HalfedgeIter h02 = h01->next();

      FaceIter f0 = h00->face();

      VertexIter v00 = h00->vertex();
      VertexIter v01 = h01->vertex();
      VertexIter v02 = h02->vertex();

    // face 1
      HalfedgeIter h10 = h00->twin();
      HalfedgeIter h11 = h10->next();
      HalfedgeIter h12 = h11->next();

      FaceIter f1 = h10->face();

      VertexIter v10 = h12->vertex();

      HalfedgeIter h01_twin = h01->twin();
      HalfedgeIter h02_twin = h02->twin();
      HalfedgeIter h11_twin = h11->twin();
      HalfedgeIter h12_twin = h12->twin();

      EdgeIter e02 = h02->edge();
      EdgeIter e11 = h11->edge();

      //check to preserve manifold
      if (v00->isBoundary() || v01->isBoundary() ||
          v00->degree() < 3 || v01->degree() < 3 ||
          e0->length() == 0.0) {
          return this->verticesEnd();
      }

      set<VertexIter> check;
      HalfedgeIter h0 = v00->halfedge();

      do{
          HalfedgeIter h_twin = h0->twin();
          check.insert(h_twin->vertex());
          h0 = h_twin->next();
      } while (h0 != v00->halfedge());

     HalfedgeIter h1 = v01->halfedge();
     int count = 0;
     do {
         HalfedgeIter h_twin = h1->twin();
         if(check.find(h_twin->vertex()) != check.end() && count == 2) {
             return this->verticesEnd();
         }
         if (check.find(h_twin->vertex()) != check.end() && count < 2) {
             count ++;
         }
         h1 = h_twin->next();
     } while (h1 != v01->halfedge());

     //time to do collapse
     //maintain all edges that connect to v00
     v00->position = (v00->position + v01->position) / 2.0;

     //halfedge changes
     //modify all edges that connect to v01
     HalfedgeIter h2 = v01->halfedge();
     do{
         EdgeIter e = h2->edge();
         HalfedgeIter h_twin = h2->twin();
         if (e != h01->edge() && e != h00->edge() && e != h12->edge()) {
             if (h2->vertex() == v01) {
                 h2->vertex() = v00;
             }
             if (h_twin->vertex() == v01) {
                 h_twin->vertex() = v00;
             }
         }
         h2 = h_twin->next();
     } while (h2 != v01->halfedge());

     //modify halfedges along e02, e11
     h01_twin->setNeighbors(h01_twin->next(), h02_twin, v02, e02, h01_twin->face());
     h02_twin->setNeighbors(h02_twin->next(), h01_twin, v00, e02, h02_twin->face());
     h11_twin->setNeighbors(h11_twin->next(), h12_twin, v10, e11, h11_twin->face());
     h12_twin->setNeighbors(h12_twin->next(), h11_twin, v00, e11, h12_twin->face());

     //vertex change
     if (v00->halfedge() == h11 || v00->halfedge() == h00) {
         v00->halfedge() = h02_twin;
     }
     if (v02->halfedge() == h02) {
         v02->halfedge() = h02_twin->next();
     }
     if (v10->halfedge() == h12) {
         v10->halfedge() = h11_twin;
     }

     //edge change
     if (e02->halfedge() == h02) {
         e02->halfedge() = h01_twin;
     }
     if(e11->halfedge() == h11) {
         e11->halfedge() = h12_twin;
     }

     //delete
     this->deleteVertex(v01);

     this->deleteEdge(h00->edge());
     this->deleteEdge(h01->edge());
     this->deleteEdge(h12->edge());

     this->deleteHalfedge(h00);
     this->deleteHalfedge(h01);
     this->deleteHalfedge(h02);
     this->deleteHalfedge(h10);
     this->deleteHalfedge(h11);
     this->deleteHalfedge(h12);

     this->deleteFace(f0);
     this->deleteFace(f1);

     return v00;
    }

  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
      
      // initialize isNew values for all existing mesh components
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->isNew = false;
      }
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          e->isNew = false;
      }
      
      // loop to set newPosition for all existing vertices
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          int n = 1;
          HalfedgeIter start = v->halfedge();
          std::vector<VertexCIter> neighbors;
          neighbors.push_back(start->next()->vertex());
          HalfedgeIter current = start->twin()->next();
          while (current != start) {
              neighbors.push_back(current->next()->vertex());
              n++;
              current = current->twin()->next();
          }
          float u;
          if (n == 3) {
              u = 3.0/16.0;
          } else {
              u = 3.0/(8.0*float(n));
          }
          Vector3D original_neighbor_position_sum = Vector3D(0.0,0.0,0.0);
          for (int i = 0; i < neighbors.size(); i++) {
              original_neighbor_position_sum += neighbors[i]->position;
          }
          v->newPosition = (1.0 - float(n) * u) * v->position + u * original_neighbor_position_sum;
      }
      
      // loop to set newPosition for future new vertices on all edges
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          HalfedgeIter h = e->halfedge();
          VertexIter A = h->vertex();
          VertexIter B = h->next()->vertex();
          VertexIter C = h->next()->next()->vertex();
          VertexIter D = h->twin()->next()->next()->vertex();
          e->newPosition = 3.0/8.0 * (A->position + B->position) + 1.0/8.0 * (C->position + D->position);
      }
      
      
      // loop to split all edges
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if ((!e->isNew) && (!e->halfedge()->vertex()->isNew && !e->halfedge()->twin()->vertex()->isNew)) {
              mesh.splitEdge(e);
          }
      }
      
      // loop through all edges to flip the appropriate ones
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if (((!e->halfedge()->vertex()->isNew && e->halfedge()->next()->vertex()->isNew) || (e->halfedge()->vertex()->isNew && !e->halfedge()->next()->vertex()->isNew)) && (e->isNew)) {
              mesh.flipEdge(e);
          }
      }
      
      // loop through all vertices to assign to new position
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
            v->position = v->newPosition;
      }

  }

  /**
   * Initializer for an edge record
   * 1. Compute a quadric for the edge as the sum of the quadrics at endpoints.
   * 2. Build a 3x3 linear system for the optimal collapsed point, as described above.
   * 3. Solve this system and store the optimal point in EdgeRecord::optimalPoint.
   * 4. Compute the corresponding error value and store it in EdgeRecord::cost.
   * 5. Store the edge in EdgeRecord::edge
   * @param edge
   * @returns an EdgeRecord
   */
  EdgeRecord::EdgeRecord( EdgeIter& _edge ) {
    // TODO: FINAL project task 3 part 1.
    Matrix4x4 K = _edge->halfedge()->vertex()->quadric + _edge->halfedge()->twin()->vertex()->quadric;
    double data[] = {K(0, 0), K(0, 1), K(0, 2),
                     K(1, 0), K(1, 1), K(1, 2),
                     K(2, 0), K(2, 1), K(2, 2)};
    Matrix3x3 A(data);
    Vector3D b(-K(0, 3), -K(1, 3), -K(2, 3));
    Vector3D x = A.inv() * b;
    optimalPoint = x;

    Vector4D u(x.x, x.y, x.z, 1);
    score = dot(u, K * u);

    edge = _edge;
  }

  Matrix4x4 calcFaceQuadric(FaceIter& f) {
    HalfedgeIter halfedge = f->halfedge();
    VertexIter A = halfedge->vertex();
    VertexIter B = halfedge->next()->vertex();
    VertexIter C = halfedge->next()->next()->vertex();
    Vector3D B_A = B->position - A->position;
    Vector3D C_A = C->position - A->position;
    Vector3D plane_normal = cross(B_A, C_A);
    double d = -dot(plane_normal, A->position); // d = -(ax + by + cz)
    Vector4D v(plane_normal.x, plane_normal.y, plane_normal.z, d);
    return outer(v, v);
  }

  Matrix4x4 calcVertexQuadric(VertexIter& v) {
    Matrix4x4 vertexQ;
    vertexQ.zero();

    HalfedgeIter h = v->halfedge();
    do {
      HalfedgeIter h_twin = h->twin();
      vertexQ += h_twin->face()->quadric;
      h = h_twin->next();
    } while(h != v->halfedge());
    return vertexQ;
  }


  /**
   * Downsampling via quadric error simplification
   * 1. Compute quadrics for each face by simply writing the plane equation for
   * that face in homogeneous coordinates, and building the corresponding quadric
   * matrix using Matrix4x4::outer(). This matrix should be stored in Face::quadric.
   *
   * 2. Compute an initial quadric for each vertex by adding up the quadrics at
   * all the faces touching that vertex. This matrix should be stored in
   * Vertex::quadric. (Note that these quadrics will get updated as edges are collapsed.)
   *
   * 3. For each edge, create an EdgeRecord and insert it into one global MutablePriorityQueue.
   *
   * 4. Until a target number of triangles is reached, collapse the best/cheapest
   * edge (as determined by the priority queue), and set the quadric at the new vertex
   * to the sum of the quadrics at the endpoints of the original edge.
   * You will also have to update the cost of any edge connected to this vertex.
   *
   * @param mesh
   */
  void MeshResampler::downsample( HalfedgeMesh& mesh )
  {
    // TODO FINAL project task 3 part 2.
    // This routine should decrease the number of triangles in the mesh using quadric error metrics.
    // STEP 1: Compute quadrics for each face
    // cout << "entering step 1" << endl;
    for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
      f->quadric = calcFaceQuadric(f);
    }
    // STEP 2: Compute an initial quadric for each vertex
    // cout << "entering step 2" << endl;
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->quadric = calcVertexQuadric(v);
    }
    // STEP 3: create an EdgeRecord and insert it into one global queue
    // cout << "entering step 3" << endl;
    MutablePriorityQueue<EdgeRecord> queue;
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      EdgeRecord eRecord(e);
      e->record = eRecord;
      queue.insert(eRecord);
    }
    // STEP 4: collapse the best/cheapest edge
    /*
     * 1. Get the cheapest edge from the queue.
     * 2. Remove the cheapest edge from the queue by calling pop().
     * 3. Compute the new quadric by summing the quadrics at its two endpoints.
     * 4. Remove any edge touching either of its endpoints from the queue.
     * 5. Collapse the edge.
     * 6. Set the quadric of the new vertex to the quadric computed in Step 2.
     * 7. Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
     */
    // cout << "entering step 4" << endl;
    Size targetNumFaces = mesh.nFaces() / 4;

    while (mesh.nVertices() > 4 and mesh.nFaces() > targetNumFaces) {
      if (queue.empty())
        break;
      /// * 1. Get the cheapest edge from the queue.
      /// * 2. Remove the cheapest edge from the queue by calling pop().
      EdgeRecord curEdgeRecord = queue.top();
      EdgeIter curEdge = curEdgeRecord.edge;
      queue.pop();

      /// 3. Compute the new quadric by summing the quadrics at its two endpoints.
      VertexIter A = curEdge->halfedge()->vertex();
      VertexIter B = curEdge->halfedge()->twin()->vertex();
      Matrix4x4 newQ = A->quadric + B->quadric;

      /// 4. Remove any edge touching either of its endpoints from the queue.
      // cout << "entering 4 A" << endl;
      HalfedgeIter h_A = A->halfedge();
      do {
        HalfedgeIter h_A_twin = h_A->twin();
        queue.remove(h_A_twin->edge()->record);
        h_A = h_A_twin->next();
      } while (h_A != A->halfedge());

      // cout << "entering 4 B" << endl;
      HalfedgeIter h_B = B->halfedge();
      do {
        HalfedgeIter h_B_twin = h_B->twin();
        queue.remove(h_B_twin->edge()->record);
        h_B = h_B_twin->next();
      } while (h_B != B->halfedge());

      queue.remove(curEdge->record);
      /// 5. Collapse the edge.
      // cout << "entering 5" << endl;
      VertexIter newVertex = mesh.collapseEdge(curEdge);
      if (newVertex == mesh.verticesEnd())
        continue;

      /// 6. Set the quadric of the new vertex to the quadric computed in Step 2.
      newVertex->quadric = newQ;
      newVertex->position = curEdgeRecord.optimalPoint;

      /// 7. Insert any edge touching the new vertex into the queue, creating new edge records for each of them.
      // cout << "entering 7" << endl;
      HalfedgeIter h = newVertex->halfedge();
      do {
        HalfedgeIter h_twin = h->twin();
        EdgeRecord newR = EdgeRecord(h_twin->edge());
        h_twin->edge()->record = newR;
        queue.insert(newR);
        h = h_twin->next();
      } while (h != newVertex->halfedge());

    }
  }

    /** Returns mean edge length of mesh for Isotropic Remeshing */
    double meanEdgeLength(HalfedgeMesh& mesh )
    {
      double sum = 0;
      double n = 0;

      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          sum += e->length();
          n += 1;
      }

      return (sum / n) * 0.9;
    }

    /** Computes the total deviation from regular degree 6 for Isotropic Remeshing */
    int deviation(int v1, int v2, int v3, int v4)
    {
        //total deviation  = |a1-6| + |a2-6| + |b1-6| + |b2-6|
        return abs(v1 - 6) + abs(v2 - 6) + abs(v3 - 6) + abs(v4 - 6);
    }

    /** Computes the average of the neighboring vertex positions and stores it in Vertex::centroid */
    void Vertex::computeCentroid() {
        //implement Vertex::computeCentroid(),
        // computes average position of the neighbors of a vertex and stores it in the member Vertex::centroid
        HalfedgeIter h = this->halfedge();
        Vector3D pos_sum = Vector3D(0, 0, 0);
        int n = 0;

        do
            {
             HalfedgeIter twin = h->twin();
             pos_sum += h->next()->vertex()->position;
             n += 1;
             h = twin->next();
             }
        while( h != this->halfedge() );

        this->centroid = pos_sum / n;
    }

    void MeshResampler::isomesh( HalfedgeMesh& mesh )
    {
        // TODO FINAL project task 4 isotropic remesh method.
        // This routine should make the mesh as "uniform" as possible
        //uniform = triangles are as close as possible to equilateral triangles of equal size, with vertex degrees as close as possible to 6

        //L = mean edge length
        double L = meanEdgeLength(mesh);
        //if an edge is longer then 4 * L / 3, split it
        double too_long = 4.0 * L / 3.0;
        //if an edge is shorter than 4 * L / 5, collapse it
        double too_short = 4.0 * L / 5.0;

        //1) If an edge is too long, split it.
        EdgeIter e = mesh.edgesBegin();
        while ( e != mesh.edgesEnd()) {
            EdgeIter next = e;
            next++;
            if (e->length() > too_long) {
                mesh.splitEdge(e);
            }
            e = next;
        }

        //2) If an edge is too short, collapse it.
        EdgeIter edge = mesh.edgesBegin();
        while (edge != mesh.edgesEnd()) {
            EdgeIter next = edge;
            next++;
            if (edge->length() < too_short) {
                EdgeIter e1 = edge->halfedge()->next()->edge();
                EdgeIter e2 = edge->halfedge()->next()->next()->edge();
                EdgeIter e3 = edge->halfedge()->twin()->next()->edge();
                EdgeIter e4 = edge->halfedge()->twin()->next()->next()->edge();
                while (next == edge || next == e1 || next == e2 || next == e3 || next == e4) {
                    next++;
                }
                mesh.collapseEdge(edge);
            }
            edge = next;
        }

        //3) If flipping an edge improves the degree of neighboring vertices, flip it.
            //want to flip an edge any time it reduces the total deviation from regular degree 6
            EdgeIter edge0 = mesh.edgesBegin();
          while (edge0 != mesh.edgesEnd()) {
              EdgeIter next = edge0;
              next++;

              HalfedgeIter h1 = edge0->halfedge();
              HalfedgeIter h2 = h1->twin();
              VertexIter v1 = h1->vertex();
              VertexIter v2 = h2->vertex();
              VertexIter v3 = h1->next()->next()->vertex();
              VertexIter v4 = h2->next()->next()->vertex();


              //total deviation in the initial configuration = |a1-6| + |a2-6| + |b1-6| + |b2-6|
              int init_dev = deviation(v1->degree(), v2->degree(), v3->degree(), v4->degree());

              //compute deviation after a "hypothetical edge flip"
              int flip_dev = deviation(v1->degree() - 1, v2->degree() - 1, v3->degree() + 1, v4->degree() + 1);

              //if the deviation decreases, then flip the edge
              if (init_dev > flip_dev) {
                  EdgeIter flipped = mesh.flipEdge(edge0);
                  if (flipped == mesh.edgesEnd()) {
                      edge0 = next;
                      continue;
                  }
                  edge0 = flipped++;
                  continue;
              }
              edge0 = next;
          }

        //4) Move vertices toward the average of their neighbors
            //mesh will have well shaped elements if each vertex is located at the center of its neighbors
            //update vertices with new positions
            //idea for new positions: rather than just replacing vertex positions by their centroids, can make more stable by moving the vertex gently toward its centroid
            VertexIter v = mesh.verticesBegin();
            while (v != mesh.verticesEnd()) {
                //p = original vertex position, c = centroid, w = weight = 1 / 5, q = p + w(c - p)
                //must move vertex only in tangent direction to avoid pushing mesh in or out, set v = v - dot(n, v)n, where n = vertex normal
                VertexIter next = v;
                next++;
                v->computeCentroid();
                Vector3D c = v->centroid;
                Vector3D n = v->normal();
                Vector3D d = c - v->position;
                d -= dot(dot(n, d), n);
                v->position += 1.0 / 5.0 * d;
                v = next;
            }
    }
}
