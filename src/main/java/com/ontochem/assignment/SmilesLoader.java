/*
 * Copyright OntoChem GmbH.
 */
package com.ontochem.assignment;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * Read smiles and ocid from input file.
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
public class SmilesLoader {
	
  private final static Pattern PAT_TAB_SPLIT = Pattern.compile( "\t" );

  private final static Logger LOG = Logger.getLogger( SmilesLoader.class.getName() );

  // ------------------------------------------------------------------------
  /**
   * Read SMILES file. 
   * 
   * @param _fileName  SMILES file; TAB separated text file with SMILES in first column,
   *                   OCID/ID in 2nd column
   * 
   * @return  map with OCID as key, SMILES as value
   * 
   * @throws IOException
   */
	public static Map<String,String> readSmiles( String _fileName ) throws IOException {
		
		LOG.info( "reading compounds: " + _fileName );
	  
		try ( BufferedReader inCsv = new BufferedReader( 
		                               new InputStreamReader(
		                                 new FileInputStream( _fileName ), "UTF8" ) ); ) {	
		  
	    final Map<String,String> targetMap = new HashMap<>();
	    
			String inLine;
			int count=0;
			while ( ( inLine = inCsv.readLine() ) != null ) {
			  
				if ( inLine.startsWith( "#" ) ) continue;
				
				String[] splitLine = PAT_TAB_SPLIT.split( inLine );
				if ( splitLine.length < 2 ) continue;
				count++;
				String smiles = splitLine[0];
				String ocid   = splitLine[1];
				
				targetMap.put( ocid, smiles  );
				
				//System.out.println( "\tread " + ocid + " " + smiles );
			}
			
			//LOG.info( "...read smiles: " + count );
	    
	    return targetMap;
		}
	}

}

