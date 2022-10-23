package com.ontochem.assignment;

import java.util.logging.Logger;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.base.exceptions.EmptyMoleculeException;
import ambit2.smarts.IsomorphismTester;
import ambit2.smarts.SmartsHelper;
import ambit2.smarts.SmartsManager;
import ambit2.smarts.SmartsParser;
import ambit2.smarts.groups.GroupMatch;
import ambit2.base.exceptions.EmptyMoleculeException;

/**
 * Author@ Shadrack Jabes., B
 * Author@ lutz.weber@ontochem.com
 * 
 * Date@ Jan 2022
 *
 * Description@
 * atom by atom search (ABAS) using cheminformatics libraries
 * 
 **/
public class StructureSearchEngine {
	
	private final static Logger LOG = Logger.getLogger( StructureSearchEngine.class.getName() );
	
	public static int searchBySubstructure( String _smiles, String _smarts, String _module, boolean _aromatic, boolean _verbose ) throws Exception {
		int rsp = 0;
		try {
			if ( _module.toLowerCase().equals( "cdk" ) ) return searchBySubstructureCdk( _smiles, _smarts );
			
			else if ( _module.toLowerCase().equals( "ambit" ) ) return searchBySubstructureAmbit( _smiles, _smarts, _aromatic, _verbose );
			
			//else if ( _module.equals( "Chemaxon" ) ) return searchBySubstructureChemaxon( _smiles, _smarts );
			
			else {
				LOG.warning( "error: chemistry module not found ");
			}
				
		} catch ( Exception e ) {
			LOG.info( "ERROR: error in searchBySubstructure module " + e );
		}
		return -1;
	}
	
	/*
	 * chemaxon substructure searcher
	 
	public static int searchBySubstructureChemaxon( String targetSmiles, String targetQuery ) {
		try  {
			String regex = "\\s+$";
			Molecule target = MolImporter.importMol(targetSmiles, "smiles");
			target.aromatize(MoleculeGraph.AROM_BASIC);
				
			Molecule query = MolImporter.importMol(targetQuery, "smarts");
			query.aromatize(MoleculeGraph.AROM_BASIC);
				
			StandardizedMolSearch ss = new StandardizedMolSearch();
			MolSearchOptions mso = ss.getSearchOptions();
		    mso.setVagueBondLevel(SearchConstants.VAGUE_BOND_OFF);
				
			ss.setTarget(target);
			ss.setQuery(query);
				
			if ( ss.isMatching() ) return 1;
			else return 0;
		    
		} catch (Exception e) {
			System.err.println( "Chemaxon error SSS: " + e) ;
			return -1;
		}
	}*/
	
	/*
	 * CDK SSS substructure searcher used to determine stereochemistry match
	 */
	public static int searchBySubstructureCdk( String _smiles, String _smarts ) throws Exception{ 
		try {
			SmilesParser smilesparser = new org.openscience.cdk.smiles.SmilesParser( SilentChemObjectBuilder.getInstance() );
			smilesparser.kekulise( true );		
			IAtomContainer mol = smilesparser.parseSmiles( _smiles );
			
			Pattern query = org.openscience.cdk.smarts.SmartsPattern.create( _smarts );
			int nUniqueHits = query.matchAll( mol ).countUnique();
			if ( nUniqueHits >0 ) return 1;
			else return 0;
			
	    } catch (Exception e) {
	    	LOG.info( "ERROR: CDK error SSS: " + _smiles + " smarts: " + _smarts );
			return -1;
		}
	}
	
	/*
	 * Ambit SSS substructure searcher
	 */
	public static int searchBySubstructureAmbit( String _smiles, String _smarts, boolean _aromatic, boolean _verbose ) { 
		try {
			
			SmartsManager man = new ambit2.smarts.SmartsManager( SilentChemObjectBuilder.getInstance() );
			IAtomContainer mol = SmilesHandler( _smiles, _aromatic ) ;  //CDK container
			 
			try {
				man.setQuery( _smarts );
				String error = man.getErrors();
				if ( error.length() > 1 ) {
					System.out.println( "Ambit smarts error: " + error );
					 System.out.println( "smarts: " + _smarts );
				}
			} catch ( Exception ee ) {
				System.out.println( "Ambit error: " + ee);
			}	
			
			if ( man.searchIn( mol ) ) {
				if ( _verbose ) System.out.println( "found: " + _smarts );
				return 1;
			} else return 0;
			
	    } catch ( Exception e ) {
	    	LOG.info( "ERROR: Ambit substructure search error: " + _smiles + " smarts: " + _smarts );
		}
		return -1;
	}
	
	/*
	 * Nick Kochev 2022-02-18 GroupMatch
	 */
	public static int searchBySubstructureAmbitAllInstances( String _smiles, String _smarts, boolean _verbose ) throws Exception {
        
		try {
        	IAtomContainer mol = SmartsHelper.getMoleculeFromSmiles( _smiles, true ) ;
        	SmartsParser sp = new SmartsParser();
        
            IsomorphismTester isoTester = new IsomorphismTester();
            isoTester.setFlagCheckStereoElements(true);
            
            GroupMatch groupMatch = new GroupMatch( _smarts, sp, isoTester );
            int posCount = groupMatch.matchCount( mol );
            
            if ( _verbose ) System.out.println( "Group " + _smarts + " found at " + posCount + " positions in " + _smiles );
            return posCount;
            
        } catch (Exception e) {
        	LOG.info( "ERROR: Ambit error all instances processing: " +_smiles + " smarts: "+_smarts);
        }
		return -1;
    }	 

	/*
	 * smiles preprocessing using CDK version 2.4.0
	 */
	public static IAtomContainer SmilesHandler( String _smiles, boolean _aromatic ) {
		try {
			IAtomContainer 	 mol 		= SmartsHelper.getMoleculeFromSmiles( _smiles, false );
			ElectronDonation model  	= ElectronDonation.daylight();
			//ElectronDonation model  	= ElectronDonation.cdk();
			CycleFinder      cycles 	= Cycles.or( Cycles.all(), Cycles.all(6) );
			Aromaticity      aroma 		= new Aromaticity( model, cycles );
			if ( _aromatic ) aroma.apply( mol );
			return  mol;
		} catch ( Exception e ) {
			LOG.info( "ERROR: Smiles Handler Error " + e );
			return null;
		}
	}
	
	/*
	 * aromatize smiles using CDK
	 */
	public static Aromaticity loadAromatizationModule() {
		try {
			ElectronDonation model = org.openscience.cdk.aromaticity.ElectronDonation.daylight();
			//ElectronDonation model = org.openscience.cdk.aromaticity.ElectronDonation.cdk();
			CycleFinder cycles = Cycles.all();
			Aromaticity arom = new org.openscience.cdk.aromaticity.Aromaticity( model, cycles );
			return arom;
		} catch (Exception e) {
			System.err.println( "load Aromatization Module Error " + e );
			return null;
		}
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

}

