import Jama.Matrix;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.StreamSupport;

import org.apache.commons.cli.*;
import org.gephi.graph.api.*;
import org.gephi.io.exporter.api.ExportController;
import org.gephi.io.exporter.spi.GraphExporter;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.EdgeDefault;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.openide.util.Lookup;
import org.gephi.graph.api.Node;

public class PMPClustering {
    public static UndirectedGraph readFromBenchmarkFile(String fileName, GraphModel graphModel)
    {
        UndirectedGraph graph = graphModel.getUndirectedGraph();
        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.split("\t");
                int vertexA = Integer.parseInt(tokens[0]);
                int vertexB = Integer.parseInt(tokens[1]);
                Node nodeA = graph.getNode(String.valueOf(vertexA));
                Node nodeB = graph.getNode(String.valueOf(vertexB));
                if ( nodeA == null)
                {
                    nodeA = graph.getGraphModel().factory().newNode(String.valueOf(vertexA));
                    nodeA.getNodeData().setLabel(String.valueOf(vertexA));
                    graph.addNode(nodeA);
                }
                if (nodeB == null)
                {
                    nodeB = graph.getGraphModel().factory().newNode(String.valueOf(vertexB));
                    nodeB.getNodeData().setLabel(String.valueOf(vertexB));
                    graph.addNode(nodeB);
                }
                graph.addEdge(nodeA, nodeB);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return graph;
    }
    /**
     * @param fileName
     * @param distanceMatrixCalculator
     * @param communitiesFile
     * @param clustersNumber
     * @param threshold
     * @param benchmark
     * @param lineGraph
     */

    public void doClustering(String fileName, DistanceMatrixCalculator distanceMatrixCalculator, EdgeClustering edgeClustering,  String communitiesFile, int clustersNumber, double threshold, boolean benchmark,  boolean lineGraph) {
        File file = new File(fileName);
        if (!file.exists()) {
            System.out.println("Error: file does not exist: " + file.toString());
            return;
        }
        String suffix = file.getName();
        suffix = suffix.split("\\.")[0];
        System.out.println("Suffix: " + suffix);
        suffix = suffix + "_" + distanceMatrixCalculator.getShortName() + "_" + edgeClustering.getShortName() + "_" + String.valueOf(clustersNumber);

        //Init a project - and therefore a workspace
        ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
        pc.newProject();
        Workspace workspace = pc.getCurrentWorkspace();

        //Get controllers and models
        ImportController importController = Lookup.getDefault().lookup(ImportController.class);
        Container container = null;
        Graph graph;
        if (benchmark) {
            graph = PMPClustering.readFromBenchmarkFile(fileName, Lookup.getDefault().lookup(GraphController.class).getModel(workspace));
            System.out.println("Edge count: " + graph.getEdgeCount());

            //for (Edge edge : graph.getEdges()) {
            //    System.out.print(edge.getSource().getNodeData().getLabel() + " ");
            //    System.out.println(edge.getTarget().getNodeData().getLabel());
            //}
        }
        else {
            //Import file
            try {
                container = importController.importFile(file.getCanonicalFile());
                container.getLoader().setEdgeDefault(EdgeDefault.UNDIRECTED);
                container.setAllowAutoNode(false);  //Don't create missing nodes
            } catch (Exception ex) {
                ex.printStackTrace();
                return;
            }

            //Append imported data to GraphAPI
            importController.process(container, new DefaultProcessor(), workspace);

            //Get a graph model - it exists because we have a workspace
            GraphModel graphModel = Lookup.getDefault().lookup(GraphController.class).getModel();
            graph = graphModel.getGraph();

            try {
                PrintWriter writer = new PrintWriter("network_format" + suffix + ".dat", "UTF-8");
                for (Edge edge : graph.getEdges()) {
                    writer.println(String.valueOf(edge.getTarget().getId()) + "\t" + String.valueOf(edge.getSource().getId()));
                }
                writer.close();
            } catch (Exception ex) {
                ex.printStackTrace();
                return;
            }

        }

        System.out.println("start distance map calculation. \n");

        Matrix d;
        d = distanceMatrixCalculator.calculateDistanceMatrix(graph);

        System.out.println("edgeDistanceMap has been calculated. \n");

        //--------------------------------------------end of matrix calculation----------------------------------------------------

        //edge clustering
        PBPolynomial sol = new PBPolynomial();

        // <Cluster, Vector of Edges>
        HashMap<Integer, Vector<Integer>> cluster_edges = edgeClustering.clusterEdges(clustersNumber, d, suffix);

        //results of clustering
        HashMap<Integer/*Node*/, HashMap<Integer, Integer>/* cluster id->count */> node_clusters = new HashMap<>();

        Vector<Integer> c = new Vector<>(cluster_edges.keySet()); //?
        c.sort(Integer::compareTo);
        System.out.println("Sorted clusters indices: " + Arrays.asList(c));
        for (Integer cluster: cluster_edges.keySet()) {
            Vector<Integer> edges = cluster_edges.get(cluster);
            for(Integer edge: edges) {
                Edge e = graph.getEdge(edge);
                Node target = e.getTarget();
                Node source = e.getSource();

                Integer targetId, sourceId;
                if(benchmark) {
                    targetId = Integer.parseInt(target.getNodeData().getLabel());
                    sourceId = Integer.parseInt(source.getNodeData().getLabel());
                }
                else {
                    targetId = target.getId();
                    sourceId = source.getId();
                }
                node_clusters.computeIfAbsent(targetId, k -> new HashMap<>());
                node_clusters.computeIfAbsent(sourceId, k -> new HashMap<>());
                if (node_clusters.get(targetId).get(c.indexOf(cluster) + 1) == null)
                    node_clusters.get(targetId).put(c.indexOf(cluster) + 1, 1);
                else
                    node_clusters.get(targetId).put(c.indexOf(cluster) + 1, node_clusters.get(targetId).get(c.indexOf(cluster) + 1) + 1 );

                if (node_clusters.get(sourceId).get(c.indexOf(cluster) + 1) == null)
                    node_clusters.get(sourceId).put(c.indexOf(cluster) + 1, 1);
                else
                    node_clusters.get(sourceId).put(c.indexOf(cluster) + 1, node_clusters.get(sourceId).get(c.indexOf(cluster) + 1) + 1 );
            }
        }

        System.out.println("PMP clusters: " + Arrays.asList(node_clusters));

        //claculate the final overlapping structure of communities
        HashMap<Integer/*Cluster id*/, Set<Integer>/*Nodes*/> format_clusters = new HashMap<>();

        for (Integer nodeId : node_clusters.keySet()) {
            if (node_clusters.get(nodeId).size() > 0) {
                double sum = 0;

                for (Integer val : node_clusters.get(nodeId).values()) {
                    sum += val;
                }
                System.out.print("Node: " + String.valueOf(nodeId) + " ");
                for (Integer clusterId : node_clusters.get(nodeId).keySet()) {
                    DecimalFormat df = new DecimalFormat("####0.00");
                    double belonging = node_clusters.get(nodeId).get(clusterId)  / sum;
                    System.out.print(String.valueOf(clusterId) + " " + df.format(belonging*100) + "% ");

                    if (belonging >= threshold) {
                        format_clusters.computeIfAbsent(clusterId,  kk -> new TreeSet<>());
                        format_clusters.get(clusterId).add(nodeId);
                    }
                }
                System.out.println();
            }
        }

        System.out.println("Clusters: " + Arrays.asList(format_clusters));

        sol.writeForMetrics(format_clusters, "pmp_" + suffix + ".dat");
        if(benchmark)
        {
            HashMap<Integer, Set<Integer>> cl = sol.formatClustersFile(communitiesFile);
            sol.writeForMetrics(cl, "truth_" + suffix + ".dat");
            System.out.println("Clusters: " + Arrays.asList(cl));
        }

        float [] r = {1.0f, 0.0f, 0.8f, 0.0f, 0, 1, 1, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f };
        float [] g = {0.0f, 0.6f, 0.0f, 0.6f, 1, 0, 1, 0.7f, 0.6f, 0.4f, 0.5f, 0.3f };
        float [] b = {0.0f, 1.0f, 0.8f, 0.2f, 1, 1, 1, 0.6f, 0.7f, 0.3f, 0.4f, 0.5f };

        Set <Node> coloredNodes = new HashSet();

        for (Node node: graph.getNodes())
            node.getNodeData().getAttributes().setValue("ClusterCount", 0);

        //gather everything in one place
        int clusterNumber = 0;
        for (Vector <Integer> cluster: cluster_edges.values()) {
            for (Integer edge_id: cluster) {
                Edge edge = graph.getEdge(edge_id);
                edge.getEdgeData().setColor(r[clusterNumber], g[clusterNumber], b[clusterNumber]);
                //edge.getEdgeData().getAttributes().setValue("Cluster", clusterNumber);

                int target_count = (int) edge.getTarget().getNodeData().getAttributes().getValue("ClusterCount");
                int source_count = (int) edge.getSource().getNodeData().getAttributes().getValue("ClusterCount");

                boolean write_target = true, write_source = true;

                if (target_count != 0) {
                    for (int j = 1; j <= target_count; j++) {
                        int clust = (int) edge.getTarget().getNodeData().getAttributes().getValue("ClusterNumb" + j);
                        if (clust == (c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1)) {
                            write_target = false;
                            break;
                        }
                    }
                }

                if (source_count != 0) {
                    for (int j = 1; j <= source_count; j++) {
                        int clust = (int) edge.getSource().getNodeData().getAttributes().getValue("ClusterNumb" + j);
                        if (clust == (c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1)) {
                            write_source= false;
                            break;
                        }
                    }
                }

                //что вот тут такое ниже - не понятно!
                edge.getEdgeData().getAttributes().setValue("ClusterNumb", (c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1));
                edge.getEdgeData().setLabel(String.valueOf(c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1));
                edge.getEdgeData().getAttributes().setValue("EdgeId", edge_id);
                if (write_source)
                {
                    source_count++;
                    edge.getSource().getNodeData().getAttributes().setValue("ClusterNumb" + source_count, (c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1));
                    edge.getSource().getNodeData().getAttributes().setValue("ClusterCount", source_count);
                }
                if (write_target)
                {
                    target_count++;
                    edge.getTarget().getNodeData().getAttributes().setValue("ClusterNumb" + target_count, (c.indexOf(cluster_edges.keySet().toArray()[clusterNumber]) + 1));
                    edge.getTarget().getNodeData().getAttributes().setValue("ClusterCount", target_count);
                }
            }
            clusterNumber++;
        }

        //Now we can assign a color to each node

        for (Node node: graph.getNodes()) {
            float rr = (float) StreamSupport.stream(graph.getEdges(node).spliterator(), false).mapToDouble(e -> e.getEdgeData().r()).sum();
            float gg = (float) StreamSupport.stream(graph.getEdges(node).spliterator(), false).mapToDouble(e -> e.getEdgeData().g()).sum();
            float bb = (float) StreamSupport.stream(graph.getEdges(node).spliterator(), false).mapToDouble(e -> e.getEdgeData().b()).sum();
            node.getNodeData().setColor(rr / graph.getDegree(node), gg / graph.getDegree(node) , bb / graph.getDegree(node));
        }

        ExportController ec = Lookup.getDefault().lookup(ExportController.class);
        GraphExporter exporter = (GraphExporter) ec.getExporter("gexf");

        if (lineGraph) {
            Workspace lineGraphWorkspace = pc.newWorkspace(pc.getCurrentProject());
            GraphModel resultLineGraphModel  = Lookup.getDefault().lookup(GraphController.class).getModel(lineGraphWorkspace);
            Graph resultLineGraph = resultLineGraphModel.getUndirectedGraph();

            HashMap <Edge, Node> edgeNodeMap = new HashMap<Edge, Node>();

            graph.getEdges().forEach(edge -> {
                Node newNode = resultLineGraphModel.factory().newNode(String.valueOf(edge.getId()));
                //покарасить в цвет изначального ребра...
                newNode.getNodeData().setLabel( edge.getEdgeData().getAttributes().getValue("ClusterNumb").toString());

                float rr = edge.getEdgeData().r();
                float gg = edge.getEdgeData().g();
                float bb = edge.getEdgeData().b();

                newNode.getNodeData().setColor(rr,gg,bb);


                resultLineGraph.addNode(newNode);
                edgeNodeMap.put(edge, newNode);
            });

            graph.getNodes().forEach(node->{
                for (Node i: graph.getNeighbors(node))
                    for (Node j: graph.getNeighbors(node))
                        if ( i.getId() >  j.getId() ) {

                            Edge orgEdgeI = graph.getEdge(node,i);
                            Edge orgEdgeJ = graph.getEdge(node,j);
                            assert orgEdgeI != null : String.format("originalEdgeI (%d,%d) is null ", node.getId(),j.getId());
                            assert orgEdgeJ != null : String.format("originalEdgeJ (%d,%d) is null ", node.getId(),j.getId());

                            Edge newEdge = resultLineGraphModel.factory().newEdge( edgeNodeMap.get(orgEdgeI), edgeNodeMap.get(orgEdgeJ) );
                            resultLineGraph.addEdge(newEdge);
                            //тут можно задать цвет ребра, например как средний
                        }
            });

            exporter.setWorkspace(lineGraphWorkspace);

            try {
                ec.exportFile(new File( String.format("%s_out_line.gexf", suffix )), exporter);
            } catch (IOException ex) {
                ex.printStackTrace();
                return;
            }

            System.out.println("The line graph has been prepared successfully!");

        } else {
            try {
                ec.exportFile(new File( String.format("%s_out.gexf", suffix) ) , exporter);
            } catch (IOException ex) {
                throw new Error(ex);
            }

        }

        System.out.println("The work has been done successfully!");
    }

    public static void main(String[] args) {
        PMPClustering clustering = new PMPClustering();
        Options options = new Options();
        Option input = new Option("i", "input", true, "path to the input file in GML format");
        input.setRequired(true);
        options.addOption(input);

        Option distanceTypeOption = new Option("d", "distance", true,
                "the type of function to measure the distance between nodes.\n" +
                        "Possible values: \nsp (shortest path) \ngd (Generalize Degree) \ncm (Commute Distance)\nacm (Amplified Commute Distance).\n" +
                        "If this option is omitted, the amplified commute distance function will be used.");
        distanceTypeOption.setRequired(false);
        options.addOption(distanceTypeOption);

        Option clusterNumberOption = new Option("k",  true, "Number of clusters to detect (int)");
        clusterNumberOption.setRequired(true);
        options.addOption(clusterNumberOption);

        Option thresholdOption = new Option("t", "threshold", true,
                "Vertex i belongs to cluster C, if node i has fraction of edges in cluster C greater then this value\n" +
                        "Default value is 0.");
        thresholdOption.setRequired(false);
        options.addOption(thresholdOption);

        Option lineGraph = new Option("l", "linegraph", false, "produce line graph as output");
        lineGraph.setRequired(false);
        options.addOption(lineGraph);

        Option edgeClusteringOption = new Option("a", "algorithm", true,
                "Specifies algorithm that will be used to find disjoint edge clusters. Possible values: \n" +
                        "pmp (p-median exact algorithm).  \n" +
                        "kmd (k-medoids heuristic)\n" +
                        "kmn (k-means heuristic).\n" +
                        "If this option is omitted, the P-Median algorithm will be used.");
        edgeClusteringOption.setRequired(false);
        options.addOption(edgeClusteringOption);

        //benchmarking
        Option benchmark = new Option("b", "benchmark", false, "use benchmark format");
        benchmark.setRequired(false);
        options.addOption(benchmark);

        //benchmarking
        Option communitiesFileOption = new Option("c", "communitiesFile", true, "use communities file to validate benchmark");
        communitiesFileOption.setRequired(false);
        options.addOption(communitiesFileOption);

        Option forceOption = new Option("f", "force", false, "" +
                "The previous founded solution by lp_solver will not be used. lp_solver will be started forcibly.\n" +
                "The flag effects only the PMP exact edge clustering algorithm.\n" +
                "By default algorithm will try to get previous founded solution");
        forceOption.setRequired(false);
        options.addOption(forceOption);

        Option output = new Option("o", "output", true,
                "the name of the output file. If option was omitted, the name of the output file\n" +
                        "will be composed automatically as \"{suffix}_{distanceName}_{clusterNumber}_out.gexf\" ");
        output.setRequired(false);
        options.addOption(output);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();

        CommandLine cmd;
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("PMPClustering <-i input_file> <-k number_of_clusters> [-ftbcd]", options);
            System.exit(1);
            return;
        }

        String inputFilePath = cmd.getOptionValue("input");
        int clusterNumber = Integer.valueOf(cmd.getOptionValue("k"));

        double threshold = 0; //default value for threshold
        if (cmd.hasOption("threshold"))
            threshold = Double.parseDouble(cmd.getOptionValue("threshold"));



        String communitiesFile = null; //"communities.dat"
        if (cmd.hasOption("benchmark") ) {
            communitiesFile = cmd.getOptionValue("communitiesFile");
        }

        DistanceMatrixCalculator distanceMatrixCalculator = new AmplifiedCommuteDistance(); //default distance calculator

        if (cmd.hasOption("distance")) {
            switch (cmd.getOptionValue("distance").toLowerCase()) {
                case "sp":
                    distanceMatrixCalculator = new ShortestPathDistance();
                    break;
                case "gd":
                    distanceMatrixCalculator = new GeneralizeDegreeDistance();
                    break;
                case "cm":
                    distanceMatrixCalculator = new CommuteDistance();
                    break;
                case "acm":
                    distanceMatrixCalculator = new AmplifiedCommuteDistance();
                    break;
                default: {
                    System.out.println("Wrong argument after distance option.");
                    System.out.println("Possible values:\n " +
                            "sp (shortest path)\n" +
                            "gd (generalize Degree)\n" +
                            "cm (commute/resistance distance)\n" +
                            "acm (amplified commute distance).");
                    System.exit(1);
                }
            }
        }

        EdgeClustering edgeClustering = new PBPolynomial();

        if (cmd.hasOption("algorithm")) {
            switch (cmd.getOptionValue("algorithm").toLowerCase()) {
                case "pmp":
                    edgeClustering = new PBPolynomial(cmd.hasOption("force"));
                    break;
                case "kmd":
                    edgeClustering = new KMedoidsWrapper();
                    break;
                case "kmn":
                    edgeClustering = new KMeansWrapper();
                    break;
                default: {
                    System.out.println("Wrong argument after edge clustering algorithm option.");
                    System.out.println("Possible values: \n" +
                            "pmp (p-median exact algorithm).  \n" +
                            "kmd (k-medoids heuristic)\n" +
                            "kmn (k-means heuristic).\n" +
                            "If this option is omitted, the P-Median algorithm will be used.");
                    System.exit(1);
                }
            }
        }

        clustering.doClustering(inputFilePath, distanceMatrixCalculator, edgeClustering, communitiesFile, clusterNumber, threshold, cmd.hasOption("benchmark"), cmd.hasOption("linegraph"));
    }
}

