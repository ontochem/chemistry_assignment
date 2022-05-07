/*
 * Copyright OntoChem GmbH.
 */
package com.ontochem.assignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Read structure ontology to perform compound assignment, creating container
 * containing is_a, has_a and querylist.
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
 * @author lutz.weber@ontochem.com
 */
public class OntologyLoader {
	
	// ==== class OntologyData ================================================
	/**
	 * Container holding parsed ontology data.
	 */
	public final static class OntologyData {
    
	    private final Map<String,Set<String>>  ocidChildMap  = new HashMap<>();
	    private final Map<String,List<String>> ocidSmartsMap = new HashMap<>();
	    private final Map<String,Set<String>>  ocidParentMap = new HashMap<>();
	    private final Map<String,String>       ocidNameMap   = new HashMap<>();
	    
	    public Map<String, Set<String>>  getOcidChildMap() { return ocidChildMap; }
	    public Map<String, String>       getOcidNameMap()  { return ocidNameMap; }
	    public Map<String, Set<String>>  getOcidParentMap() { return ocidParentMap; }
	    public Map<String, List<String>> getOcidSmartsMap() { return ocidSmartsMap; }
	    
	    public void setOcidChildren( String _id, Set<String> _children ) {
	    	ocidChildMap.put( _id, _children );
	    }
	    public void setOcidName( String _id, String _name ) {
	    	ocidNameMap.put( _id, _name );
	    }
	    public void setOcidParents( String _id, Set<String> _parents ) {
	    	ocidParentMap.put( _id, _parents );
	    }
	    public void setOcidSmarts( String _id, List<String> _smarts ) {
	    	ocidSmartsMap.put( _id, _smarts );
	    }
	}
	// ========================================================================
	  
		
	private final static Logger LOG = Logger.getLogger( OntologyLoader.class.getName() );
  
	/**
	 * Reads/parses compound classes OBO with SMARTS.
	 * 
	 * @param _inObo
	 * @param _module  lower case name of chemical library used (e.g. cdk, ambit, chemaxon)
	 * 
	 * @return
	 * 
	 * @throws IOException
	 */
	public static OntologyData readObo( String _inObo, String _module ) throws IOException {
		
		LOG.info( "reading OBO: " + _inObo );
		
		final OntologyData ontData = new OntologyData();

		final boolean isModuleCdkOrAmbit = ChemLib.CHEMLIB_CDK.equals( _module ) || ChemLib.CHEMLIB_AMBIT.equals( _module );
		final boolean isModuleCA         = ChemLib.CHEMLIB_CA.equals( _module );

		try ( BufferedReader inObo = new BufferedReader( 
		                               new InputStreamReader( 
		                                 new FileInputStream( new File( _inObo ) ), "UTF8" ) ); ) {
		
	  		String inLine = null;
	  		while ( ( inLine = inObo.readLine() ) != null ) {
	  			
	  			if ( inLine.startsWith( "[Term]" ) ) {
			      
	  				final List<String> smartsList	= new ArrayList<>();
	  				final Set<String>  childSet	= new HashSet<>();
	  				final Set<String>  parentSet	= new HashSet<>();
	  				String name = null;
	  				String id   = null;
	
	  				String inConceptLine = null;
	  				while ( ( inConceptLine = inObo.readLine() ) != null ) {
	
				        if ( inConceptLine.startsWith( "[" ) ) {
				        	throw new IOException( "Start of new Stanza without previous empty line." );
				        }
			        
				        if ( inConceptLine.trim().isEmpty() ) {
				        	break;
				        }
			        
				        int tagSepOff = inConceptLine.indexOf( ':' );
				        if ( tagSepOff < 2 ) {
				        	continue;
				        }
					    final String tag   = inConceptLine.substring( 0, tagSepOff );
					    String value = inConceptLine.substring( tagSepOff + 1 ).trim();
			        
					    int commentSepOff = value.indexOf( " !" );
					    if ( commentSepOff >= 0 ) {
					    	value = value.substring( 0, commentSepOff ).trim();
					    }
			        
				        if ( "id".equals( tag ) ) {
				        	id = value;
				        } else if ( "name".equals( tag ) ) {
				        	name = value;
					    } else if ( "is_a".equals( tag ) ) {
					    	parentSet.add( value ); 
				        } else if ( "has_a".equals( tag ) ) {
				        	childSet.add( value );
				        } else if ( tag.endsWith( "_smarts" ) ) {
				        	if ( ( isModuleCdkOrAmbit && "cdk_smarts".equals( tag ) ) ||
				               ( isModuleCA && "oc_smarts".equals( tag ) ) ) {
				        		String smarts = value.replace("\\!","!").replace("\\\\","\\");
				        		smartsList.add( smarts );
				        	}
				        }
	  				}
			      
	  				if ( id != null ) {
			            //leaves without smarts are useless!!!
			            if ( ! ( ( ! parentSet.isEmpty() ) && smartsList.isEmpty() && childSet.isEmpty() ) ) {
			  		        if ( name != null ) ontData.setOcidName( id, name );
			  		        ontData.setOcidChildren( id, childSet );
			  		        ontData.setOcidParents( id, parentSet );
			  		        ontData.setOcidSmarts( id, smartsList );
			            }
	  				}
	  			}
	  		}
		}
		return ontData;
	}
	
	// ==== command line usage (testing) ======================================
	public static void main( String[] _args ) throws Exception {
		readObo( _args[0], _args[1] );
	}
}
