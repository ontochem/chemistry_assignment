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

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Read smiles and ocid from input file.
 * 
 * <h3>Changelog</h3>
 * <ul>
 *   <li>2022-05-15
 *     <ul>
 *       <li>third version</li>
 *     </ul>
 *   </li>
 * </ul>
 * 
 * @author lutz.weber@ontochem.com
 */
public class SmilesLoader {
	
	private final static Pattern PAT_TAB_SPLIT = Pattern.compile( "\t" );

	private final static Logger LOG = Logger.getLogger( SmilesLoader.class.getName() );

	/**
	* Read SMILES file. 
	* 
	* @param _fileName  SMILES file; TAB separated text file with SMILES in first column,
	*                   OCID/ID in 2nd column
	* 
	* @return  map with OCID as key, SMILES as value
	 * @throws Exception 
	*/
	public static Map<String,String> readSmiles( String _fileName, int _max, boolean _verbose ) throws Exception {
		
		LOG.info( "reading compounds: " + _fileName );
		final Map<String,String> targetMap = new HashMap<>();
		
		SmilesGenerator _smilesG = getCdkSmilesGenerator();
		SmilesParser _smilesP = getCdkSmilesParser();
		Aromaticity _arom = loadAromatizationModule();
		
		try ( BufferedReader inCsv = new BufferedReader( 
		                               new InputStreamReader(
		                                 new FileInputStream( _fileName ), "UTF8" ) );){	
			String inLine;
			int count=0;
			
			while ( ( inLine = inCsv.readLine() ) != null ) {
				if ( inLine.startsWith( "#" ) ) continue;
				count++;
				if ( count > _max ) continue;
				String[] splitLine = inLine.split( "\t" );
				if ( splitLine.length < 2 ) continue;
				String smiles = splitLine[0];
				String ocid   = splitLine[1];
				try {
					String normSmiles = normalizeSmilesCDK( smiles, _smilesG, _smilesP, _arom );
					if ( _verbose ) System.out.println ( normSmiles );
					targetMap.put( ocid, normSmiles  );
				} catch (Exception e) {
					// TODO: handle exception
				}
			}
			LOG.info( "...read smiles: " + count );
		} catch ( Exception er ) {
			LOG.info( "ERROR: " + er.getLocalizedMessage() );
		}
		return targetMap;
		
	}
	
	public static String normalizeSmilesCDK( String _smiles, SmilesGenerator _smilesG, 
												SmilesParser _smilesP, Aromaticity _arom ) throws Exception {
		
		IAtomContainer mol = _smilesP.parseSmiles( _smiles );
		_arom.apply( mol );
		org.openscience.cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens( mol ); 
		String stdSmiles = _smilesG.create( mol );
		
        return stdSmiles;
    }
	
	public static org.openscience.cdk.smiles.SmilesGenerator getCdkSmilesGenerator() {
		org.openscience.cdk.smiles.SmilesGenerator smilesG = 
			new org.openscience.cdk.smiles.SmilesGenerator( SmiFlavor.Isomeric );
		return smilesG;
	}
	
	public static org.openscience.cdk.smiles.SmilesParser getCdkSmilesParser() {
		org.openscience.cdk.smiles.SmilesParser smilesP = 
			new org.openscience.cdk.smiles.SmilesParser( SilentChemObjectBuilder.getInstance() );
		smilesP.kekulise( true );
		return smilesP;
	}
	
	public static Aromaticity loadAromatizationModule() {
		try {
			ElectronDonation model = org.openscience.cdk.aromaticity.ElectronDonation.daylight();
			//ElectronDonation model = org.openscience.cdk.aromaticity.ElectronDonation.cdk();
			CycleFinder cycles = Cycles.all();
			Aromaticity arom = new org.openscience.cdk.aromaticity.Aromaticity( model, cycles );
			return arom;
		} catch (Exception e) {
			LOG.info( "ERROR: load Aromatization Module Error " + e );
			return null;
		}
	}

}

