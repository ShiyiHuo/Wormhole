import edu.princeton.cs.algs4.SeparateChainingHashST;
import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.EdgeWeightedDigraph;
import edu.princeton.cs.algs4.DirectedEdge;
import edu.princeton.cs.algs4.EdgeWeightedDirectedCycle;
import edu.princeton.cs.algs4.Stack;
import java.util.Iterator;
import java.util.NoSuchElementException;


public class Worm {
	private int V;   // # of vertices in the graph
    	public int W;   // # of wormholes in the graph
	private DirectedEdge[] wormholeArray;
	private SeparateChainingHashST<String,Vertex> hash1;
	private SeparateChainingHashST<Integer,Vertex> hash2;
	public AdjMatrixEdgeWeightedDigraph G;
	public FloydWarshall f;
    

	public Worm( In in ) { 
		// Create a new problem instance.
		V=in.readInt();
		hash1=new SeparateChainingHashST<String,Vertex>(V);
		hash2=new SeparateChainingHashST<Integer,Vertex>(V);
		G=new AdjMatrixEdgeWeightedDigraph(V);
		
		//build hash table
		for (int i=0; i<V; i++) {
			String name=in.readString();
			int x_coor=in.readInt();
			int y_coor=in.readInt();
			int z_coor=in.readInt();
			Vertex v=new Vertex(i,name,x_coor,y_coor,z_coor);
			hash1.put(name,v);
			hash2.put(i,v);
		}
		
        //read wormholes and set up an array of wormholes
        W=in.readInt();
        wormholeArray=new DirectedEdge[W];
        for (int i=0; i<W; i++) {
            String entrance=in.readString();
            int startPoint=hash1.get(entrance).getIndex();
            String exit=in.readString();
            int endPoint=hash1.get(exit).getIndex();
            wormholeArray[i]=new DirectedEdge(startPoint,endPoint,0);
        }
        
		//add edges to graph
		for (int s=0; s<V; s++) {
			for (int t=0; t<V; t++) {
				if (s!=t) {
                    if (W!=0) {
                        for (int k=0; k<W; k++) {
                            //if there is a wormhole between two vertices
                            if (wormholeArray[k].from()==s && wormholeArray[k].to()==t) {
                                DirectedEdge e=new DirectedEdge(s,t,0);
                                G.addEdge(e);
                            }
                        }
                    }
                    double weight=calculateDistance(hash2.get(s).getName(),hash2.get(t).getName());
                    DirectedEdge e=new DirectedEdge(s,t,weight);
                    G.addEdge(e);	
				}
			}
		}
        
		//run Floyd on graph G
		f=new FloydWarshall(G); 
	}

	public double dist( String origP, String destP ) {
		// return the distance from origP to destP
		return f.dist(hash1.get(origP).getIndex(), hash1.get(destP).getIndex());
	}
	
	//calculate Euclidean distance between two vertices
	public double calculateDistance (String origP, String destP) {
		double dist=0.0;
		Vertex v=hash1.get(origP);
		int v_xcoor=v.getx_coor();
		int v_ycoor=v.gety_coor();
		int v_zcoor=v.getz_coor();
		Vertex w=hash1.get(destP);
		int w_xcoor=w.getx_coor();
		int w_ycoor=w.gety_coor();
		int w_zcoor=w.getz_coor();
		dist=Math.sqrt(Math.pow(w_xcoor-v_xcoor,2)+Math.pow(w_ycoor-v_ycoor,2)+Math.pow(w_zcoor-v_zcoor,2));
		return dist;
	}
    
    
    // least number of wormholes in any shortest path from origP to destP
    // Get the path of floyd and check how many edges on the path has weight 0,
    // increment wormNum if an edge weighed 0 is found.
	public int worms( String origP, String destP ) {
        int[][] wormArrayFloyd=f.wormArray();
        return wormArrayFloyd[hash1.get(origP).getIndex()][hash1.get(destP).getIndex()];
	}
 
	public String query( String origP, String destP ) {
		// Output the "The distance from ... wormholes." string.
		return "The distance from "+origP+" to "+destP+" is "+Math.round(dist(origP,destP))+" using "+worms(origP,destP)+" worm holes.";
	}
  
	public static void main(String[] args) {
		// You can test your program with something like this.
		In in = new In( args[0] );
		int T = in.readInt();
		for (int t=1; t<=T; t++) {
			System.out.println("Case " + t + ":") ;
			Worm w = new Worm( in );
			int Q = in.readInt() ;
			for (int i=0; i<Q; i++) {
				String p1s = in.readString() ;
				String p2s = in.readString() ;
				System.out.println( w.query( p1s, p2s ) ) ;
			}
		}
		
	}
}


class Vertex {
	private int index;
	private String name;
	private int x_coor;
	private int y_coor;
	private int z_coor;
	
	public Vertex (int index, String name, int x_coor, int y_coor, int z_coor) {
		this.index=index;
		this.name=name;
		this.x_coor=x_coor;
		this.y_coor=y_coor;
		this.z_coor=z_coor;
	}
    
    public int getIndex() {
        return this.index;
    }
    
    public String getName() {
        return this.name;
    }
    
    public int getx_coor() {
        return this.x_coor;
    }
    
    public int gety_coor() {
        return this.y_coor;
    }
    
    public int getz_coor() {
        return this.z_coor;
    }
	
}


class FloydWarshall {
    private boolean hasNegativeCycle;  // is there a negative cycle?
    private double[][] distTo;  // distTo[v][w] = length of shortest v->w path
    private DirectedEdge[][] edgeTo;  // edgeTo[v][w] = last edge on shortest v->w path
    private int[][] worms;  //needs to initialize it with worms in the graph
    private DirectedEdge[][] graphMatrix;
    
    /**
     * Computes a shortest paths tree from each vertex to to every other vertex in
     * the edge-weighted digraph <tt>G</tt>. If no such shortest path exists for
     * some pair of vertices, it computes a negative cycle.
     * @param G the edge-weighted digraph
     */
    public FloydWarshall(AdjMatrixEdgeWeightedDigraph G) {
        int V = G.V();
        distTo = new double[V][V];
        edgeTo = new DirectedEdge[V][V];
        worms=new int[V][V];
        graphMatrix=G.adjMatrix();
        
        // initialize distances to infinity and worms array according to weight of edges
        for (int v = 0; v < V; v++) {
            for (int w = 0; w < V; w++) {
                distTo[v][w] = Double.POSITIVE_INFINITY;
                if (v!=w && graphMatrix[v][w].weight()==0.0) {
                    worms[v][w]=1;
                }
            }
        }

        // initialize distances using edge-weighted digraph's
        for (int v = 0; v < G.V(); v++) {
            for (DirectedEdge e : G.adj(v)) {
                distTo[e.from()][e.to()] = e.weight();
                edgeTo[e.from()][e.to()] = e;
            }
            // in case of self-loops
            if (distTo[v][v] >= 0.0) {
                distTo[v][v] = 0.0;
                edgeTo[v][v] = null;
            }
        }
        
        // Floyd-Warshall updates
        for (int i = 0; i < V; i++) {
            // compute shortest paths using only 0, 1, ..., i as intermediate vertices
            for (int v = 0; v < V; v++) {
                if (edgeTo[v][i] == null) continue;  // optimization
                for (int w = 0; w < V; w++) {
                    if (distTo[v][w] > distTo[v][i] + distTo[i][w]) {
                        distTo[v][w] = distTo[v][i] + distTo[i][w];
                        edgeTo[v][w] = edgeTo[i][w];
                        worms[v][w] = worms[v][i] + worms[i][w];
                    }
                    
                    if (distTo[v][w] == distTo[v][i] + distTo[i][w]) {
                        if (worms[v][w]>worms[v][i]+worms[i][w]) {
                            edgeTo[v][w]=edgeTo[i][w];
                            worms[v][w]=worms[v][i]+worms[i][w];
                        }
                    }
                }
                // check for negative cycle
                if (distTo[v][v] < 0.0) {
                    hasNegativeCycle = true;
                    return;
                }
            }
        }
    }
    
    public int[][] wormArray() {
        return worms;
    }
    
    public boolean hasNegativeCycle() {
        return hasNegativeCycle;
    }
    
    public Iterable<DirectedEdge> negativeCycle() {
        for (int v = 0; v < distTo.length; v++) {
            // negative cycle in v's predecessor graph
            if (distTo[v][v] < 0.0) {
                int V = edgeTo.length;
                EdgeWeightedDigraph spt = new EdgeWeightedDigraph(V);
                for (int w = 0; w < V; w++)
                    if (edgeTo[v][w] != null)
                        spt.addEdge(edgeTo[v][w]);
                EdgeWeightedDirectedCycle finder = new EdgeWeightedDirectedCycle(spt);
                assert finder.hasCycle();
                return finder.cycle();
            }
        }
        return null;
    }
    
    /**
     * Is there a path from the vertex <tt>s</tt> to vertex <tt>t</tt>?
     * @param s the source vertex
     * @param t the destination vertex
     * @return <tt>true</tt> if there is a path from vertex <tt>s</tt>
     * to vertex <tt>t</tt>, and <tt>false</tt> otherwise
     */
    public boolean hasPath(int s, int t) {
        return distTo[s][t] < Double.POSITIVE_INFINITY;
    }
    
    /**
     * Returns the length of a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>.
     * @param s the source vertex
     * @param t the destination vertex
     * @return the length of a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>;
     * <tt>Double.POSITIVE_INFINITY</tt> if no such path
     * @throws UnsupportedOperationException if there is a negative cost cycle
     */
    public double dist(int s, int t) {
        if (hasNegativeCycle())
            throw new UnsupportedOperationException("Negative cost cycle exists");
        return distTo[s][t];
    }
    
    /**
     * Returns a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>.
     * @param s the source vertex
     * @param t the destination vertex
     * @return a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>
     * as an iterable of edges, and <tt>null</tt> if no such path
     * @throws UnsupportedOperationException if there is a negative cost cycle
     */
    public Iterable<DirectedEdge> path(int s, int t) {
        if (hasNegativeCycle())
            throw new UnsupportedOperationException("Negative cost cycle exists");
        if (!hasPath(s, t)) return null;
        Stack<DirectedEdge> path = new Stack<DirectedEdge>();
        for (DirectedEdge e = edgeTo[s][t]; e != null; e = edgeTo[s][e.from()]) {
            path.push(e);
        }
        return path;
    }
    
    // check optimality conditions
    private boolean check(EdgeWeightedDigraph G, int s) {
        
        // no negative cycle
        if (!hasNegativeCycle()) {
            for (int v = 0; v < G.V(); v++) {
                for (DirectedEdge e : G.adj(v)) {
                    int w = e.to();
                    for (int i = 0; i < G.V(); i++) {
                        if (distTo[i][w] > distTo[i][v] + e.weight()) {
                            System.err.println("edge " + e + " is eligible");
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }
    
}


class AdjMatrixEdgeWeightedDigraph {
    private static final String NEWLINE = System.getProperty("line.separator");
    
    private int V;
    private int E;
    private DirectedEdge[][] adj;

    public AdjMatrixEdgeWeightedDigraph(int V) {
        if (V < 0) throw new RuntimeException("Number of vertices must be nonnegative");
        this.V = V;
        this.E = 0;
        this.adj = new DirectedEdge[V][V];
    }

    public DirectedEdge[][] adjMatrix() {
        return adj;
    }
    
    public int V() {
        return V;
    }

    public int E() {
        return E;
    }
    
    public void addEdge(DirectedEdge e) {
        int v = e.from();
        int w = e.to();
        if (adj[v][w] == null) {
            E++;
            adj[v][w] = e;
        }
    }
    
    public Iterable<DirectedEdge> adj(int v) {
        return new AdjIterator(v);
    }
    
    // support iteration over graph vertices
    private class AdjIterator implements Iterator<DirectedEdge>, Iterable<DirectedEdge> {
        private int v;
        private int w = 0;
        
        public AdjIterator(int v) {
            this.v = v;
        }
        
        public Iterator<DirectedEdge> iterator() {
            return this;
        }
        
        public boolean hasNext() {
            while (w < V) {
                if (adj[v][w] != null) return true;
                w++;
            }
            return false;
        }
        
        public DirectedEdge next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            return adj[v][w++];
        }
        
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(V + " " + E + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (DirectedEdge e : adj(v)) {
                s.append(e + "  ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }

}
