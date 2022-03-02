package com.ontochem.assignment;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
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
	
	public static int searchBySubstructure( String _smiles, String _smarts, String module) throws Exception {
		int rsp = 0;
		try {
			if ( module.equals( "Cdk" ) ) return searchBySubstructureCdk( _smiles, _smarts );
			
			if ( module.equals( "Ambit" ) ) return searchBySubstructureAmbit( _smiles, _smarts );
			
			if ( module.equals( "Chemaxon" ) ) return searchBySubstructureChemaxon( _smiles, _smarts );
				
		} catch ( Exception e ) {
			System.out.println( "error searchBySubstructure: " + e );
		}
		return -1;
	}
	
	/*
	 * chemaxon substructure searcher
	 */
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
	}
	
	/*
	 * CDK SSS substructure searcher
	 */
	public static int searchBySubstructureCdk( String _smiles, String _smarts ) throws Exception{ 
		try {
			SmilesParser smilesparser = new org.openscience.cdk.smiles.SmilesParser( SilentChemObjectBuilder.getInstance() );
			smilesparser.kekulise( true );			
			IAtomContainer target = smilesparser.parseSmiles( _smiles );
			org.openscience.cdk.smarts.SmartsPattern query = org.openscience.cdk.smarts.SmartsPattern.create( _smarts );
			int nUniqueHits = query.matchAll( target ).countUnique();
			return nUniqueHits;
			
	    } catch (Exception e) {
			System.out.println( "CDK error SSS: " + _smiles + " smarts: " + _smarts );
			return -1;
		}
	}
	
	/*
	 * Ambit SSS substructure searcher
	 */
	public static int searchBySubstructureAmbit( String _smiles, String _smarts ) { 
		try {
			IAtomContainer mol = SmilesHandler( _smiles ) ;  //CDK container
			SmartsManager man = new ambit2.smarts.SmartsManager( SilentChemObjectBuilder.getInstance() );
			man.setUseCDKIsomorphismTester( false );
			try {
				man.setQuery( _smarts );
				String error = man.getErrors();
				if ( error.length() > 1 ) {
					System.out.print( "Ambit smarts error: " + error );
					System.out.println( "smarts: " + _smarts );
				}
			} catch ( Exception ee ) {
				System.out.println( "Ambit error: " + ee);
			}	
			if ( man.searchIn( mol ) ) return 1;
			else return 0;
	    } catch ( Exception e ) {
			System.out.println( "Ambit substructure search error: " + _smiles + " smarts: " + _smarts );
		}
		return -1;
	}
	
	/*
	 * Nick Kochev 2022-02-18 GroupMatch
	 */
	public static int searchBySubstructureAmbitAllInstances( String _smiles, String _smarts ) throws Exception {
        
		try {
        	IAtomContainer mol = SmartsHelper.getMoleculeFromSmiles(_smiles, true) ;
        	SmartsParser sp = new SmartsParser();
            IsomorphismTester isoTester = new IsomorphismTester();
            GroupMatch groupMatch;
			
            groupMatch = new GroupMatch( _smarts, sp, isoTester );
            int posCount = groupMatch.matchCount(mol);
            
            //System.out.println( "Group " + _smarts + " found at " + posCount + " positions in " + _smiles );
            return posCount;
            
        } catch (Exception e) {
            System.out.println( "Ambit error all instances processing: " +_smiles + " smarts: "+_smarts);
        }
		return -1;
    }	 

	/*
	 * smiles preprocessing using CDK version 2.4.0
	 */
	public static IAtomContainer SmilesHandler( String _smiles ) {
		try {
			IAtomContainer 	 mol 	= SmartsHelper.getMoleculeFromSmiles( _smiles, false );
			ElectronDonation model  = ElectronDonation.daylight();
			CycleFinder      cycles = Cycles.or( Cycles.all(), Cycles.all(6) );
			Aromaticity      aroma 	= new Aromaticity( model, cycles );
			aroma.apply( mol );
			AtomContainerManipulator.convertImplicitToExplicitHydrogens( mol );
			return  mol;
		} catch ( Exception e ) {
			System.err.println( "Smiles Handler Error " + e );
			return null;
		}
	}
	
	/*
	 * aromatize smiles using CDK
	 */
	public static Aromaticity loadAromatizationModule() {
		try {
			ElectronDonation model = org.openscience.cdk.aromaticity.ElectronDonation.daylight();
			CycleFinder cycles = Cycles.all();
			Aromaticity arom = new org.openscience.cdk.aromaticity.Aromaticity(model,cycles);
			return arom;
		} catch (Exception e) {
			System.err.println( "load Aromatization Module Error " + e );
			return null;
		}
	}

}

