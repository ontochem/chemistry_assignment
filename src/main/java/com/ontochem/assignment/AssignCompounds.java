/*
* Copyright OntoChem GmbH, 
* Blücherstrasse 24, 06120 Halle (Saale), Germany
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*    http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
package com.ontochem.assignment;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.ontochem.assignment.OntologyLoader.OntologyData;

/**
 * Compound assignment to chemical structure classes in a ontology
 * version 1.06
 *
 * <h3>Changelog</h3>
 * <ul>
 *   <li>2022-10-23
 *     <ul>
 *       <li>sixth version</li>
 *       <li>includes stereochemistry assignement via CDK</li>
 *       <li>includes smarts assignement via Ambit</li>
 *     </ul>
 *   </li>
 * </ul>
 * 
 * @author shadrack.j.barnabas@ontochem.com
 * @author lutz.weber@ontochem.com
 */
public class AssignCompounds {  
	
	private final static Logger LOG = Logger.getLogger( AssignCompounds.class.getName() );
	
	/**
	 * Container for assignment parameters.
	 */
	public static class AssignmentParameters {
		
		/** lower case name of chemical library used, e.g. ambit, cdk, chemaxon,
	     *  see {@link ChemLib} */
	    private String  module;
	    private String  ontologyFilename;
	    private String  smilesFilename;
	    private String  outFilename;
	    private String  statisticsFilename;
	    private int     nThreads = 1;
	    private int 	max = 1;    //maximum number of evaluated smiles, for testing
	    private boolean appendModuleInfoToFilename = false;
	    private boolean writeToStandardOut = false;
	    private boolean writeLeavesOnly = true;
    
	    // ---- getter/setter -----------------------------------------
	    public String getModule() { return module; }
	    public AssignmentParameters setModule( String _module ) {
	    	module = trimWithEmptyAsNull( _module );
	    	if ( module != null ) module = ChemLib.resolveChemLib( module );
	    	return this;
	    }
    
	    public int getnThreads() { return nThreads; }
	    public AssignmentParameters setnThreads( int _nThreads ) {
	    	if ( _nThreads > 0 ) {
	    		nThreads = _nThreads;
	    	}
	    	return this;
	    }
	    
	    public String getOntologyFilename() { return ontologyFilename; }
	    public AssignmentParameters setOntologyFilename( String _ontologyFilename ) {
	    	ontologyFilename = trimWithEmptyAsNull( _ontologyFilename ); return this;
	    }
    
	    public boolean isAppendModuleInfoToFilename() { return appendModuleInfoToFilename; }
	    public AssignmentParameters setAppendModuleInfoToFilename( boolean _appendModuleInfoToFilename ) {
	    	appendModuleInfoToFilename = _appendModuleInfoToFilename; return this;
	    }
    
	    public String getOutFilename() { return outFilename; }
	    public AssignmentParameters setOutFilename( String _outFilename ) {
	    	outFilename = trimWithEmptyAsNull( _outFilename ); return this;
	    }
	    
	    public String getSmilesFilename() { return smilesFilename; }
	    public AssignmentParameters setSmilesFilename( String _smilesFilename ) {
	    	smilesFilename = trimWithEmptyAsNull( _smilesFilename ); return this;
	    }
	    
	    public int getMax() { return max; }
	    public AssignmentParameters setMax( int _max ) {
	    	max = _max; return this;
	    }
    
	    public String getStatisticsFilename() { return statisticsFilename; }
	    public AssignmentParameters setStatisticsFilename( String _statisticsFilename ) {
	    	statisticsFilename = trimWithEmptyAsNull( _statisticsFilename ); return this;
	    }
    
	    public boolean isWriteToStandardOut() { return writeToStandardOut; }
	    public AssignmentParameters setWriteToStandardOut( boolean _writeToStandardOut ) {
	    	writeToStandardOut = _writeToStandardOut; return this;
	    }
	    
	    public boolean isWriteLeavesOnly() { return writeLeavesOnly; }
	    public AssignmentParameters setWriteLeavesOnly( boolean _writeLeavesOnly ) {
	    	writeLeavesOnly = _writeLeavesOnly; return this;
	    }
    
	    /**
	     * Throws an {@link IOException} if any required parameter is not set.
	     * 
	     * @throws IOException
	     */
	    public void checkParameters() throws IOException {
	    	checkNotNull( "module", module );
	    	checkNotNull( "ontologyFilename", ontologyFilename );
	    	checkNotNull( "outFilename", outFilename );
	    	checkNotNull( "smilesFilename", smilesFilename );
	    }
    
	    @Override
	    public String toString() {
	    	try {
	    		return new ObjectMapper().writeValueAsString( this );
	    	} catch ( IOException ioe ) {
	    		return "ERROR serializing AssignmentParameters to JSON: " + ioe.getMessage(); 
	    	}
	    }
    
	    private void checkNotNull( String _parameter, String _value ) throws IOException {
	    	if ( _value == null ) {
	    		throw new IOException( "Parameter '" + _parameter + "' not set" );
	    	}
	    }
    
	    private String trimWithEmptyAsNull( String _str ) {
	    	String str = _str;
	    	if ( str != null ) {
	    		str = str.trim();
	    		if ( str.isEmpty() ) {
	    			str = null;
	    		}
	    	}
	    	return str;
	    }
	}
  
	/**
	 * Execute compound assignment.
	 * 
	 * @param _parameters
	 * 
	 * @throws Exception
	 */
	public static void runAssignment( AssignmentParameters _parameters ) throws Exception {
		
		_parameters.checkParameters();
		
		LOG.info( "using chemistry module: " + _parameters.getModule() );
		
		long startTime = System.nanoTime();
		boolean aromatic = true;
		boolean verbose = false;

	    /*
	     * step 1: read chemistry ontology with smarts
	     */
	    OntologyData ontData = OntologyLoader.readObo( _parameters.getOntologyFilename(), _parameters.getModule(), aromatic );
	    final Map<String,Set<String>>  ocidClass2childMap   		= ontData.getOcidChildMap();
	    final Map<String,Set<String>>  ocidClass2parentMap  		= ontData.getOcidParentMap();
	    final Map<String,List<String>> ocidClass2smartsList			= ontData.getOcidSmartsMap();
	    final Map<String,String>       ocidClass2nameMap  			= ontData.getOcidNameMap();
	    final Map<String,Set<String>>  ocidClass2offspringsMap 	    = new HashMap<String,Set<String>>();
	    final Map<String,Set<String>>  ocid2ancestorsMap 			= new HashMap<String,Set<String>>();
	    
	    if ( ocidClass2childMap.size() == 0 ) {
			LOG.info("generating ocidClass2childMap as has_a relationship was not found in the OBO ...");
			List<String> idList1 = new ArrayList<>();
	   		ocidClass2parentMap.forEach( ( key, value )->{ idList1.addAll( value ); } );
			for ( String id : idList1 ) {
				List<String> list1 = new ArrayList<>();
				Iterator<String> itre = ocidClass2parentMap.keySet().iterator();
				while ( itre.hasNext() ) {
					String key = itre.next();
					Set<String> parents = ocidClass2parentMap.get(key);
					if( parents.contains( id ) ) list1.add( key );
				}
				HashSet<String> hsList1 = new HashSet<String>( list1 );
				ocidClass2childMap.put( id, hsList1 );
			}
		}
	    
	    /* 
		 * step 2: identify root class id 
		 */
	   	String rootId = null;
	    int countRoots = 0;
	    for ( String ocid : ocidClass2parentMap.keySet() ) {
	    	if ( ocidClass2parentMap.get( ocid ).isEmpty() ) {
	    		rootId = ocid;
	    		countRoots++;
	    	}
	    }
	    if ( countRoots > 1 ) {
	    	throw new IOException( "error - more than 1 root concept found in compound classes OBO: " + countRoots );
	    } else if ( rootId == null ) {
	    	throw new IOException( "no root concept found in compound classes OBO" );
	    } 
	    
	    /*
	     * step 3: read and load smiles
	     */
	    final Map<String,String> toProcessOcid2SmilesMap = SmilesLoader.readSmiles( _parameters.getSmilesFilename(), _parameters.getMax(), verbose );
	    final List<String>       toProcessOcidList       = new ArrayList<String>( toProcessOcid2SmilesMap.keySet() );
	    
	    /*
	     * step 4: calculate ancestorMap and offspringMap
	     */
	    final Map<String,Set<String>>  ocidClass2AllAncestorsMap  = new HashMap<String, Set<String>>();
	    for ( String ocidClass : ocidClass2nameMap.keySet() ) {
	    	Set<String> ancestorSet = ancestors( ocidClass, ocidClass2parentMap );
	    	ocidClass2AllAncestorsMap.put( ocidClass, ancestorSet );
	    }
	    for ( String ocidClass : ocidClass2nameMap.keySet() ) {
	    	Set<String> offspringSet = offsprings( ocidClass, ocidClass2childMap );
	    	ocidClass2offspringsMap.put( ocidClass, offsprings( ocidClass, ocidClass2childMap ) );
	    }
	    
	    /* 
	     * step 5: follow hierarchy of class id top down and check assignment, write into ocidAssignmentMap
	     */
	    final Map<String,Set<String>> ocidAssignmentMap = 
	    		AssignmentUtils.hierarchicalParallelClassAssignment( rootId, _parameters.getModule(), aromatic, verbose,
                                                              		 _parameters.getnThreads(), 
                                                              		 toProcessOcidList, 
                                                              		 toProcessOcid2SmilesMap, 
                                                              		 ocidClass2smartsList, 
                                                              		 ocidClass2parentMap,
                                                              		 ocidClass2childMap 	);

	    String outfile = _parameters.getOutFilename();
	    if ( _parameters.isAppendModuleInfoToFilename() ) outfile = outfile + "_" + _parameters.getModule()+".tsv";
	    
	    try ( Writer out = new OutputStreamWriter(
                             new BufferedOutputStream(
                               new FileOutputStream( new File( outfile ) ) 
                             ),
                           StandardCharsets.UTF_8 ); ) {

	    	Set<String> ocidClassSet1 = new HashSet<String>();
	    	Set<String> ocidClassSet2 = new HashSet<String>();
      
		    for ( String ocid : ocidAssignmentMap.keySet() ) {
		    	
		    	ocidClassSet1.clear();
		        ocidClassSet2.clear();
		        
		        final Set<String> ocidClassSet = ocidAssignmentMap.get( ocid );
		        
		        // check ancestors, omit concept if a parent is missing
		        for ( String ocidClass : ocidClassSet ) {
		        	
		        	final Set<String> classAllAncestorsSet = ocidClass2AllAncestorsMap.get( ocidClass );
  
		        	// leave out if a concept has no parents
		        	if ( ocidClass2parentMap.get( ocidClass ).isEmpty() ) continue;
  
		        	// leave out if a ancestor concept is missing
		        	//final boolean missingParent = classAllAncestorsSet.stream()
                    //                         .anyMatch( ancestor -> ! ocidClassSet.contains( ancestor ) );
		        	boolean missingParent = false;
					for ( String ancestor : classAllAncestorsSet ) {
						if ( !ocidClassSet.contains( ancestor ) ) {
							//System.out.println( ocidClass+" missing parent: "+ancestor);
							missingParent = true;
							continue;
						}
	    			}
					
		        	if ( missingParent || ( ocidClass2smartsList.get( ocidClass ).isEmpty() ) ) {
		        		continue;
		        	} else {
		        		ocidClassSet1.add( ocidClass );
		        	}
		        }
		        
		        // check for children present, omit concept if valid child with smarts has been found
		        for ( String ocidClass : ocidClassSet1 ) {
		        	//System.out.println( ocid + " assigned to 1: "+ocidClass );
		        	boolean validChild = false;
		        	final Set<String> classAllOffspringsSet = ocidClass2offspringsMap.get( ocidClass );
		  
		        	// leave out if offspring with smarts is present, one is enough
		        	if ( classAllOffspringsSet != null ) {
		        		for ( String offspring : classAllOffspringsSet ) {
		        			if ( ocidClassSet1.contains( offspring ) && ocidClass2smartsList.containsKey( offspring ) ) {
		        				validChild = true;
		        				break;
		        			}
		        		}
		        	}
		        	if ( validChild ) {
		        		continue;
		        	} else {
		        		ocidClassSet2.add( ocidClass );
		        	}
		        }
		        
		        if ( _parameters.isWriteToStandardOut() ) {
		        	System.out.print( toProcessOcid2SmilesMap.get( ocid )+ "\t" + ocid );
		        }
		        
		        out.append( ocid + "\t" + toProcessOcid2SmilesMap.get( ocid ) +"\n");
		        
		        if ( _parameters.isWriteLeavesOnly() ) {
			        for ( String newClass : ocidClassSet2 ) {
			        	out.append( "is_a\t" ).append( newClass ).append( "\t" )
	                            .append( ocidClass2nameMap.get( newClass ) ).append( "\n" );
			        	out.flush();
			        	if ( _parameters.isWriteToStandardOut() ) {
			        		System.out.print( "\t" + newClass + "\t" + ocidClass2nameMap.get( newClass ) );
			        	}
			        }
			        if ( _parameters.isWriteToStandardOut() ) System.out.print( "\n" );
		        	
		        } else {
		        	for ( String newClass : ocidClassSet1 ) {
			        	out.append( "is_a\t" ).append( newClass ).append( "\t" )
	                            .append( ocidClass2nameMap.get( newClass ) ).append( "\n" );
			        	out.flush();
			        	if ( _parameters.isWriteToStandardOut() ) {
			        		System.out.print( "\t" + newClass + "\t" + ocidClass2nameMap.get( newClass ) );
			        	}
			        }
		        }
		        System.out.println();
		        out.write( "\n" );
		        System.out.print("\n" );
		        ocid2ancestorsMap.put( ocid, ocidClassSet1 );
		        ocidClassSet1 = new HashSet<String>();
		        ocidClassSet2 = new HashSet<String>();
		    }  
		    
		    if ( _parameters.getStatisticsFilename() != null ) {
		    	LOG.info( "statistics file name: " + _parameters.getStatisticsFilename() );
		    	
		    	Map<String,Integer> statisticsMap = new HashMap<String,Integer>();
		    	
		    	for ( String ocid : ocid2ancestorsMap.keySet() ) {
		    		Set<String> classIdSet = ocid2ancestorsMap.get( ocid );
		    		for ( String classId : classIdSet ) {
		    			if ( statisticsMap.get( classId ) != null ) {
		    				int count = statisticsMap.get( classId ) + 1;
		    				statisticsMap.put( classId, count );
		    			} else statisticsMap.put( classId, 1 );
		    		}
		    	}
		    	
		    	try ( Writer outStat = new OutputStreamWriter(
                    new BufferedOutputStream(
                      new FileOutputStream( new File( _parameters.getStatisticsFilename() ) ) 
                    ), StandardCharsets.UTF_8 ); ) {
			    		for ( String classId : statisticsMap.keySet() ) {
				        	outStat.append( classId ).append( "\t" ).append( ocidClass2nameMap.get( classId )+"\t").append( String.valueOf( statisticsMap.get( classId ) )).append( "\n" );
				        }
		    	}
	    	
		    } else LOG.info( "no statistics file written.");
	    }
	    long duration = System.nanoTime() - startTime;
    
	    LOG.info( "Elapsed time (s): " + TimeUnit.NANOSECONDS.toSeconds( duration ) );
	}
    
	// ------------------------------------------------------------------------
	
	/**
	 * Get all ancestor IDs for provided compound class concept id.
	 * 
	 * @param _classId
	 * @param _parentMap
	 * 
	 * @return
	 */
	private static Set<String> ancestors( String _classId, Map<String,Set<String>> _parentMap ) {
		
		final Set<String> ancestorList = new HashSet<>();
		
		Set<String> idList = Collections.singleton( _classId );
		
		while ( ! idList.isEmpty() ) {
	    final Set<String> newIdList = new HashSet<>();  //list of all ids for one hierarchy level
			for ( String id : idList ) {
			  final Set<String> parentList = _parentMap.get( id );
				if ( ( parentList != null ) && ( ! parentList.isEmpty() ) ) {
				  newIdList.addAll( parentList );
				}
			}
			idList = newIdList;
			ancestorList.addAll( newIdList ); //list of all ids for down of classId
		}
		return ancestorList;
	}
	
	// ------------------------------------------------------------------------
	/**
	 * Returns all descendants of provided compound class concept.
	 * 
	 * @param _classId
	 * @param _childMap
	 * 
	 * @return
	 */
	private static Set<String> offsprings( String _classId, Map<String,Set<String>> _childMap ){
		
    final Set<String> offspringList = new HashSet<>();
    
		Set<String> idList = Collections.singleton( _classId );
		
		while ( ! idList.isEmpty() ) {
			final Set<String> newIdList = new HashSet<>(); 	//list of all ids for one hierarchy level
			for ( String id : idList ) {
			  Set<String> childList = _childMap.get( id );
				if ( childList != null ) {
				  newIdList.addAll( childList );
				}
			}
			idList = newIdList;
			offspringList.addAll( newIdList ); //list of all ids for down of classId
		}
		
		return offspringList;
	}
	
	// ==== command line usage ================================================
	// ------------------------------------------------------------------------
	
	/**
	 * Prints command line usage help and exits.
	 * 
	 * @param _exitCode
	 */
	private static void usage( int _exitCode ) {
    
		System.err.println( "Run compound assignment.\n" +
                        "usage: java " + AssignCompounds.class.getName() + " OPTIONS\n" +
                        "   -t    --threads      THREADCOUNT\n" +
                        "                          number of threads to use for assignment\n" +
                        "   -max                 MAX maximal number of smiles to process\n" +
                        "   -m    --module       CHEMLIB\n" +
                        "                          chemical library to use, one of\n" +
                        "                            Cdk\n" +
                        "                            Ambit\n" +
                        "                            ChemAxon\n" +
                        "   -c    --obo          FILENAME\n" +
                        "                          name of OBO file with chemical classes\n" +
                        "   -s    --smiles-id    FILENAME\n" +
                        "                          name of file with smiles and ids\n" +
                        "   -o   --output-file   FILENAME\n" +
                        "                          name of output file the assignments will be written to\n" + 
                        "   -a   --ancestors     ANCESTORS\n" +
                        "                          all - creates output with all assigned classes\n" + 
                        "                          parents - creates output with assigned parent classes\n" + 
                        "   -sta  --statistics   FILENAME\n" +
                        "                          creates output with all assigned compounds to classes\n" + 
                        "   -ter  --terminal       write also to standard out\n" +
                        "   -u    --append-module  append module name to output-file\n" +
                        "   -h    --help           print this help and exit\n" 
                      );
    
		System.exit( _exitCode );
	}
  
	/**
	 * Parses command line parameters. In case of an error it will be written
	 * to std-err and {@link #usage(int)} is called which exits with code 1.
	 * 
	 * @param _args
	 * 
	 * @return  parsed parameters
	 */
	private static AssignmentParameters parseCommandLine( String[] _args ) {
    
		final AssignmentParameters parameters = new AssignmentParameters();
    
		for ( int argIdx = 0; argIdx < _args.length; argIdx++ ) {
      
			final String arg     = _args[ argIdx ];
			argIdx++;
			
			final String nextArg = argIdx < _args.length ? _args[ argIdx ] : null;
      
			if ( "-c".equals( arg ) || "--obo".equals( arg ) ) {
				parameters.setOntologyFilename( nextArg );
			} else if ( "-t".equals( arg ) || "--threads".equals( arg ) ) {
		        try {
		        	parameters.setnThreads( Integer.parseInt( nextArg ) );
		        } catch ( NumberFormatException nfe ) {
		        	System.err.println( "Thread count is not an integer: '" + nextArg + "'" );
		        	usage( 1 );
		        }
			} else if ( "-m".equals( arg ) || "--module".equals( arg ) ) {
				parameters.setModule( nextArg );
			} else if ( "-s".equals( arg ) || "--smiles-id".equals( arg ) ) {
				parameters.setSmilesFilename( nextArg );
			} else if ( "-o".equals( arg ) || "--output-file".equals( arg ) ) {
				parameters.setOutFilename( nextArg );
			} else if ( "-a".equals( arg ) || "--ancestors".equals( arg ) ) {
				if ( nextArg.equals("all") ) parameters.setWriteLeavesOnly(false);
				else if ( nextArg.equals("parents") ) parameters.setWriteLeavesOnly(true);
				else {
					System.err.println( nextArg+" Unknown parameter '" + arg + "'" );
					usage( 1 );
				}
			} else if ( "-sta".equals( arg ) || "--statistics".equals( arg ) ) {
				parameters.setStatisticsFilename( nextArg );
			} else if ( "-ter".equals( arg ) || "--terminal".equals( arg ) ) {
				parameters.setWriteToStandardOut( true );
				argIdx = argIdx - 1;
			} else if ( "-max".equals( arg ) || "--maximal".equals( arg ) ) {
				parameters.setMax( Integer.parseInt( nextArg ) );
			} else if ( "-app".equals( arg ) || "--append-module".equals( arg ) ) {
				parameters.setAppendModuleInfoToFilename( true );
				argIdx = argIdx - 1;
			} else if ( "-h".equals( arg ) || "--help".equals( arg ) ) {
				usage( 0 );
				argIdx = argIdx - 1;
			} else {
				System.err.println( nextArg + " Unknown parameter '" + arg + "'" );
				usage( 1 );
			}
		}
		return parameters;
	}
  
  	/**
  	 * Command line entry point.
  	 * 
  	 * @param _args
  	 * 
  	 * @throws Exception
  	 */
	public static void main( String[] _args ) throws Exception {

		AssignmentParameters parameters = parseCommandLine( _args );
    
		runAssignment( parameters );
	}

}

