/*
 * Copyright OntoChem GmbH.
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.logging.Logger;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.ontochem.assignment.OntologyLoader.OntologyData;



/**
 * Compound assignment to chemical structure classes in a ontology
 *
 * <h3>Changelog</h3>
 * <ul>
 *   <li>2022-02-25
 *     <ul>
 *       <li>initial version</li>
 *     </ul>
 *   </li>
 * </ul>
 * 
 * @author Shadrack Jabes., B
 * @author lutz.weber@ontochem.com
 */
public class AssignCompounds {  
	
  // ==== class AssignmentParameters ========================================
  /**
   * Container for assignment parameters.
   */
  public static class AssignmentParameters {
    
    private String  ontologyFilename;
    /** lower case name of chemical library used, e.g. ambit, cdk, chemaxon,
     *  see {@link ChemLib} */
    private String  module;
    private String  smilesFilename;
    private String  outFilename;
    private int     nThreads = 1;
    private boolean appendModuleInfoToFilename = false;
    private boolean writeToStandardOut = false;
    
    // ---- getter/setter -----------------------------------------
    public String getModule() { return module; }
    public AssignmentParameters setModule( String _module ) {
      module = trimWithEmptyAsNull( _module );
      if ( module != null ) {
        module = ChemLib.resolveChemLib( module );
      }
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
    
    public String getOutFilename() { return outFilename; }
    public AssignmentParameters setOutFilename( String _outFilename ) {
      outFilename = trimWithEmptyAsNull( _outFilename ); return this;
    }
    
    public String getSmilesFilename() { return smilesFilename; }
    public AssignmentParameters setSmilesFilename( String _smilesFilename ) {
      smilesFilename = trimWithEmptyAsNull( _smilesFilename ); return this;
    }
    
    public boolean isAppendModuleInfoToFilename() { return appendModuleInfoToFilename; }
    public AssignmentParameters setAppendModuleInfoToFilename( boolean _appendModuleInfoToFilename ) {
      appendModuleInfoToFilename = _appendModuleInfoToFilename; return this;
    }
    
    public boolean isWriteToStandardOut() { return writeToStandardOut; }
    public AssignmentParameters setWriteToStandardOut( boolean _writeToStandardOut ) {
      writeToStandardOut = _writeToStandardOut; return this;
    }
    
    // ------------------------------------------------------------
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
    
    // ------------------------------------------------------------
    @Override
    public String toString() {
      try {
        return new ObjectMapper().writeValueAsString( this );
      } catch ( IOException ioe ) {
        return "ERROR serializing AssignmentParameters to JSON: " + ioe.getMessage(); 
      }
    }
    
    // ------------------------------------------------------------
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
  // ========================================================================
  
  
  private final static Logger LOG = Logger.getLogger( AssignCompounds.class.getName() );
  
	/**
	 * Execute compound assignment.
	 * 
	 * @param _parameters
	 * 
	 * @throws Exception
	 */
  public static void runAssignment( AssignmentParameters _parameters ) throws Exception {

    _parameters.checkParameters();

    long startTime = System.nanoTime();

    /*
     * step 1: read and load smarts and chemistry ontology
     */
    OntologyData ontData = OntologyLoader.readObo( _parameters.getOntologyFilename(),
                                                   _parameters.getModule() );

    final Map<String,Set<String>>  ocidClass2ChildMap   = ontData.getOcidChildMap();
    final Map<String,Set<String>>  ocidClass2ParentMap  = ontData.getOcidParentMap();
    final Map<String,List<String>> ocidClass2SmartsList = ontData.getOcidSmartsMap();
    final Map<String,String>       ocidClass2NameMap  	= ontData.getOcidNameMap();
    final Map<String,Set<String>>  ocidClass2AllOffspringsMap = new HashMap<>();

    /*
     * step 2: read and load smiles
     */
    final Map<String,String> toProcessOcid2SmilesMap = SmilesLoader.readSmiles( _parameters.getSmilesFilename() );
    final List<String>       toProcessOcidList       = new ArrayList<>( toProcessOcid2SmilesMap.keySet() );

    /* 
     * step 3: identify root class id 
     */
    String rootId = ocidClass2ChildMap.keySet()
                    .stream()
                    .filter( key -> ! ocidClass2ParentMap.containsKey( key ) )
                    .findFirst().orElse( null );

    if ( rootId == null ) {
      throw new IOException( "No root concept found in compound classes OBO" );
    }

    /*
     * step 4: calculate ancestorMap and offspringMap
     */
    final Map<String,Set<String>>  ocidClass2AllAncestorsMap  = new HashMap<>();
    
    for ( String ocidClass : ocidClass2NameMap.keySet() ) {
      ocidClass2AllAncestorsMap.put( ocidClass, ancestors( ocidClass, ocidClass2ParentMap ) );
    }
    
    for ( String ocidClass : ocidClass2NameMap.keySet() ) {
      ocidClass2AllOffspringsMap.put( ocidClass, offsprings( ocidClass, ocidClass2ChildMap ) );
    }

    /* 
     * step 5: follow hierarchy of class id top down and check assignment, write into ocidAssignmentMap
     */
    final Map<String,Set<String>>
    ocidAssignmentMap = AssignmentUtils
                        .hierarchicalParallelClassAssignment( rootId, _parameters.getModule(),
                                                              _parameters.getnThreads(), 
                                                              toProcessOcidList, 
                                                              toProcessOcid2SmilesMap, 
                                                              ocidClass2SmartsList, 
                                                              ocidClass2ChildMap );

    try ( Writer out = new OutputStreamWriter(
                          new BufferedOutputStream(
                              new FileOutputStream( new File( _parameters.isAppendModuleInfoToFilename() ?
                                  _parameters.getOutFilename() + 
                                  "_" + _parameters.getModule() + ".tsv" :
                                    _parameters.getOutFilename() ) ) 
                              ),
                          StandardCharsets.UTF_8 ); ) {

      final Set<String> ocidClassSet1 = new HashSet<>();
      final Set<String> ocidClassSet2 = new HashSet<>();
      
      for ( String ocid : ocidAssignmentMap.keySet() ) {
  
        ocidClassSet1.clear();
        ocidClassSet2.clear();
  
        final Set<String> ocidClassSet = ocidAssignmentMap.get( ocid );
  
        // check ancestors, omit concept if a parent is missing
        for ( String ocidClass : ocidClassSet ) {
  
          final Set<String> classAllAncestorsSet = ocidClass2AllAncestorsMap.get( ocidClass );
  
          // leave out if a concept has no parents
          if ( ! ocidClass2ParentMap.containsKey( ocidClass ) ) continue;
  
          // leave out if a ancestor concept is missing
          final boolean 
          missingParent = classAllAncestorsSet.stream()
                                              .anyMatch( ancestor -> ! ocidClassSet.contains( ancestor ) );
          
          if ( missingParent || ( ! ocidClass2SmartsList.containsKey( ocidClass ) ) ) {
            continue;
          } else {
            ocidClassSet1.add( ocidClass );
          }
        }
  
        // check for children present, omit concept if valid child with smarts has been found
        for ( String ocidClass : ocidClassSet1 ) {
          
          boolean validChild = false;
          final Set<String> classAllOffspringsSet = ocidClass2AllOffspringsMap.get( ocidClass );
  
          // leave out if offspring with smarts is present, one is enough
          if ( classAllOffspringsSet != null ) {
            for ( String offspring : classAllOffspringsSet ) {
              if ( ocidClassSet1.contains( offspring ) && ocidClass2SmartsList.containsKey( offspring ) ) {
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
          System.out.println( ocid + "\t" + toProcessOcid2SmilesMap.get( ocid ) );
        }
        for ( String newClass : ocidClassSet2 ) {
          out.append( "\t" ).append( newClass ).append( "\t" )
                            .append( ocidClass2NameMap.get( newClass ) ).append( "\n" );
          if ( _parameters.isWriteToStandardOut() ) {
            System.out.println( ocid + "\t" + newClass + "\t" + ocidClass2NameMap.get( newClass ) );
          }
        }
        out.write( "\n" );
      }  // for
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
                        "    -t  --threads  THREADCOUNT\n" +
                        "                     number of threads to use for assignment\n" +
                        "    -m  --module  CHEMLIB\n" +
                        "                     chemical library to use, one of\n" +
                        "                       Cdk\n" +
                        "                       Ambit\n" +
                        "                       ChemAxon\n" +
                        "    -c  --obo  FILENAME\n" +
                        "                     name of OBO file with chemical classes\n" +
                        "    -s  --smiles-id  FILENAME\n" +
                        "                     name of file with smiles and ids\n" +
                        "    -o  --output-file  FILENAME\n" +
                        "                     name of output file the assignments will be written to\n" + 
                        "    -h  --help\n" +
                        "                     print this help and exit\n" 
                      );
    
    System.exit( _exitCode );
  }
  
  // ------------------------------------------------------------------------
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
      final String nextArg = argIdx + 1 < _args.length ? _args[ argIdx + 1 ] : null;
      
      
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
        
      } else if ( "-h".equals( arg ) || "--help".equals( arg ) ) {
        usage( 0 );
        
      } else {
        System.err.println( "Unknown parameter '" + arg + "'" );
        usage( 1 );
      }
    }
    
    return parameters;
  }
  
  // ------------------------------------------------------------------------
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

