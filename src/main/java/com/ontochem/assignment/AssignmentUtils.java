package com.ontochem.assignment;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.logging.Logger;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;

/**
 * Utilities for performing compound assignment.
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
public class AssignmentUtils {
	
	private final static Logger LOG = Logger.getLogger( AssignmentUtils.class.getName() );
  
	/**
	 * @param rootClassId
	 * @param _module
	 * @param _nThreads
	 * @param _ocidListIn
	 * @param _ocid2smilesMap
	 * @param _ocidClass2smartsList
	 * @param _ocidClass2parentMap
	 * @param _ocidClass2childMap
	 * 
	 * @return
	 * 
	 * @throws IOException
	 */
	public final static Map<String,Set<String>> 
	hierarchicalParallelClassAssignment( String rootClassId, String _module, boolean _aromatic, boolean _verbose,
											int _nThreads, 
											List<String> _ocidListIn, 
											Map<String,String> _ocid2smilesMap, 
											Map<String,List<String>> _ocidClass2smartsList, 
											Map<String,Set<String>> _ocidClass2parentMap,
											Map<String,Set<String>> _ocidClass2childMap ) throws IOException {
    	
		final Map<String,Set<String>> ocid2classMap = new ConcurrentHashMap<>();
		final ForkJoinPool            forkJoinPool  = new ForkJoinPool( _nThreads );
		
		try {
			
			forkJoinPool.submit( () -> {
			  
				_ocidListIn.parallelStream().forEach( ocid -> {
					
					String smiles = _ocid2smilesMap.get( ocid );
					Set<String> classIdList = new HashSet();
					classIdList.add( rootClassId );
					
					Set<String> childrenSet = new HashSet<String>();
					Set<String> assignedSet = new HashSet<>();
					
					int count    = 0;
					int countAss = 0;
					int countAll = 0;
					
					while ( !classIdList.isEmpty() ) {
					  
						final Set<String> newIdList = new HashSet<String>(); 	//list of all ids for one hierarchy level
						
						for ( String classId : classIdList ) {
							
							if ( assignedSet.contains( classId ) ) continue;
							
							boolean parentNotAssigned = false;
							Set<String> classParents = _ocidClass2parentMap.get( classId );
							if ( classParents != null ) {
								for ( String parent : classParents ) {
									if ( parentNotAssigned ) break;
									if ( !assignedSet.contains( parent ) ) parentNotAssigned = true;
								}
							}
							if ( parentNotAssigned ) continue;
							
							List<String> smartsList = _ocidClass2smartsList.get( classId );
							if ( ( smartsList != null ) && !smartsList.isEmpty() ) {
								count++;
								if ( assign( smiles, smartsList, _module, _aromatic, _verbose ) ) {
									countAss++;
									assignedSet.add( classId );
									countAll++;
								} 
							} else {
								if ( !_ocidClass2childMap.isEmpty() ) assignedSet.add( classId );
								countAll++;
							}
							
							childrenSet = _ocidClass2childMap.get( classId );
							if ( ( childrenSet != null ) && !childrenSet.isEmpty() ){
								newIdList.addAll( childrenSet );
							}
						}
						classIdList = newIdList;
					}
					ocid2classMap.put( ocid, assignedSet );
					if ( _verbose ) System.out.println( ocid + " " + assignedSet.size() );
				});
			}).get();
			
			forkJoinPool.shutdown();
		
		} catch ( Exception e ) {
			throw new IOException( "Error in parallel hierarchical class assignment: " + e.getMessage(), e );
		}
		
		return ocid2classMap;
  }
    
	/**
	 * @param _smiles
	 * @param _smartsList
	 * @param _module  lower case name of chemical library used (e.g. cdk, ambit,
	 *                 chemaxon, see {@link ChemLib})
	 * 
	 * @return
	 */
	public static boolean assign( String _smiles, List<String> _smartsList, String _module, 
															boolean _aromatic, boolean _verbose ) {
		  
		if ( ChemLib.CHEMLIB_CA.equals( _module ) ) {
			LOG.info("ChemAxon is not implemented in public version...stopping");
			return false;
			//return assignChemaxon( _smiles, _smartsList, _module );
		  
		} else if ( ChemLib.CHEMLIB_CDK.equals( _module ) || 
		            ChemLib.CHEMLIB_AMBIT.equals( _module ) ) {
		  
			return assignCdkOrAmbit( _smiles, _smartsList, _module, _aromatic, _verbose );
		  
		} else {
			LOG.severe( "Unexpected chemical library module: '" + _module + "'" );
			return false;
		}
	}
	
	/**
	 * handle logical operators in smarts and uses atom by atom search (ABAS) to assign a compound. 
	 */
	public static boolean assignCdkOrAmbit( String _smiles, List<String> _smartsList, String _module, 
																boolean _aromatic, boolean _verbose ) {
		try {
			List<String> ORResponse = new ArrayList<String>();
			List<String> NOTResponse = new ArrayList<String>();
			
			boolean isORstructure 	= false;
			boolean isNOTstructure 	= false;
			String  trueResponse 	= new String( "true" );
			int size = _smartsList.size();
		
			/*
			 * A class consists of smarts set.
			 * Smarts set is a collection of smarts query or smarts structures with logical operators
			 * the following additional logic has been added to for CDK smarts:
			 * 
			 * AND structures: query1XXXquery2, query1.query2. Both dots or XXX separated smarts are taken apart and we check each of them individually and only when all are ok then the whole expression is considered as AND connected structure.
			 * POLY structures: 3ZZZquery2 query 2 must occor 3 times
			 * NOT structures: !query, !query1XXXquery2, !query1.query2
			 * OR structures: queries without NOT symol, and AND structures after resolving XXX, WWW, dots.
			 */
			for ( int i=0; i < _smartsList.size(); i++ ) {
				
				String Qry = _smartsList.get( i );
				String Qry2 = Qry.replaceAll( "^!", "" );
				
	    		if ( Qry2.equals( Qry ) ) {								//OR smarts
	    			isORstructure = true;
	    			
	    			String Qry1 = Qry.replaceAll("XXX", "");
	    			
	    			if ( Qry1.equals( Qry ) ) {							//OR smarts
	    				
	    				String Qry3 = Qry.replace(".", "");
		    			
	    				if ( Qry3.equals( Qry ) ) {						//OR smarts, potential EXACT and MORE
			    			int rspd = checkQueryMultiplicity( _smiles, Qry, _module, _aromatic, _verbose );
			    			if ( rspd >0 ) ORResponse.add( "true" ); else ORResponse.add("false");
		    				//System.out.println("checkQueryMultiplicity: " + " OR " + rspd + " " + Qry);
			    			
	    				} else {										//OR + ANDdot smarts, potential EXACT and MORE										
 		    				List<String> newQryList = Segmenter( Qry, "." );
		    				List<String> AndConnectedList = new ArrayList<String>();
		    				for ( int j = 0; j<newQryList.size(); j++) { 
		    					String newQry = newQryList.get(j);
		    					int rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
		    					//System.out.println("sss: " + rspd + " " + newQry);
		    					if ( rspd >0 ) AndConnectedList.add("true"); else AndConnectedList.add("false");
		    				}
		    				if ( IsAndConnected( AndConnectedList ) ) { 
		    					String rspdS = AndConnectedList.get( 0 );
		    					ORResponse.add( rspdS );
		    					//System.out.println("AND with dots " + rspd + " " + Qry);
		    				}
		    			}
	    				
		    		} else {											//OR + AND smarts
		    			List<String> newQryList = Segmenter( Qry, "XXX" );
	    				List<String> AndConnectedList= new ArrayList<String>();
	    				
	    				String newQry = newQryList.get(0);				//stereospecific query preceeds before XXX
	    				int rspd = checkQueryStereoSpecificity( _smiles, newQry, _module, _aromatic, _verbose );
	    				if ( rspd >0 ) AndConnectedList.add( "true" ); else AndConnectedList.add( "false" );
	    				//System.out.println("stereo: " + rspd + " " + newQry );
	    				
	    				newQry = newQryList.get( 1 );					//non-stereospecific query follows after XXX
	    				////rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
	    				rspd = checkQueryMultiplicity( _smiles, newQry, _module, _aromatic, _verbose );
	    				if ( rspd >0 ) AndConnectedList.add( "true" ); else AndConnectedList.add( "false" );
	    				//System.out.println("sss: " + rspd + " " + newQry );
	    				
	    				if ( IsAndConnected( AndConnectedList ) ) { 
	    					String rspdS = AndConnectedList.get(0);		//only true if both part1 and part2 are true
	    					ORResponse.add( rspdS );
	    					//System.out.println("AND with XXX " + rspd + " " + Qry);
	    				}
	    			}
	    			
	    		} else {												//NOT smarts
	    			
	    			isNOTstructure = true;
	    			
	    			String Qry1 = Qry2.replaceAll("XXX", "");
	    			
	    			if( Qry1.equals( Qry2 ) ) {							//NOT smarts
	    				
	    				String Qry3 = Qry2.replace(".", "");
		    			
	    				if( Qry3.equals( Qry2 ) ) {						//NOT as simple smarts
	    					
		    				// no patterns found, OR structure preceeded by NOT
		    				int rspd1 = checkQueryMultiplicity( _smiles, Qry2, _module, _aromatic, _verbose );
		    				if ( rspd1 >0 ) NOTResponse.add("false");
	    					else NOTResponse.add("true");
			    			
		    			} else {										//NOT+AND smarts
		    				
 		    				List<String> newQryList = Segmenter( Qry2,"." );
		    				List<String> AndConnectedList= new ArrayList<String>();
		    				
		    				for (int j = 0; j<newQryList.size(); j++) {
		    					String newQry = newQryList.get(j);
		    					int rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
		    					//System.out.println("sss: " + rspd + " " + newQry);
		    					if ( rspd >0 ) AndConnectedList.add("true");
		    					else AndConnectedList.add("false");
		    				}
		    				if ( IsAndConnected( AndConnectedList ) ) { 
		    					String rspd1 = AndConnectedList.get(0);
		    					if( rspd1.equals("true") ) NOTResponse.add("false"); else NOTResponse.add("true");
		    				}
		    			}
	    				
	    			} else { 											//NOT + AND smarts
	    				List<String> newQryList = Segmenter( Qry2, "XXX" );
	    				List<String> AndConnectedList= new ArrayList<String>();
	    				
	    				for ( int j = 0; j<newQryList.size(); j++ ) {
	    					String newQry = newQryList.get(j);
	    					int rspd = 0;
	    					if ( j==0 ) {
		    					//stereospecific query preceeds before XXX
	    						rspd = checkQueryStereoSpecificity( _smiles, newQry, _module, _aromatic, _verbose );
	    					} else {
		    					//non-stereospecific query follows after XXX
	    						rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
	    					}
	    					if ( rspd >0 ) AndConnectedList.add( "true" ); else AndConnectedList.add( "false" );
	    				}
	    				if ( IsAndConnected(AndConnectedList ) ) { 
	    					String rspd1 = AndConnectedList.get(0);
	    					if ( rspd1.equals( "true" ) ) NOTResponse.add( "false" ); else NOTResponse.add( "true" );
	    				}
	    			}
	    		}
			}
			
			/*
			 * a substructure match is considered, if at least one OR structure is a match in the smarts set
			 * and all NOT structures in the smarts set did not match
			 */
			boolean asgn1 = IsOrConnected( ORResponse, "true") ;
			
			/*
			 * two scenarios in NOT structures: 1. only NOT structure 2. both NOT and OR structures
			 * scenario 1: only NOT
			 */
			if( isNOTstructure ) {
				boolean asgn2 = false;
				boolean check = IsAndConnected( NOTResponse );
				if ( check ) asgn2 = Boolean.valueOf ( NOTResponse.get(0) );
				
				if ( !isORstructure ) {
					if( asgn2 ) return true;
					else return false;
				}
				
				/*
				 * two scenarios in NOT structures: 1. only NOT structure 2. both NOT and OR structures
				 * scenario 2: both NOT and OR
				 */
				if ( isORstructure ) {
		    		if ( asgn2 ) {
		    			if ( asgn1 ) return true;
		    			else return false;
		    		}
		    	}
			}
			
			/*
			 * scenario 3: no NOT structures, only OR structures
			 */
			else if ( isORstructure ) {
		    	//System.out.println("found OR structures");
		    	if ( asgn1 ) return true;
		    	else return false;
		    }
				
		} catch (Exception e) {
			System.out.println("assigner subroutine error" + e );
		}
		return false;
	}
		
	/**
	 * @param smi
	 * @param sma
	 * @param module
	 * @return true or false. 
	 */     
    public static int checkQueryStereoSpecificity( String _smi, String _sma, String _module, 
    															boolean _aromatic, boolean _verbose ) {
    	// stereochemistry is implemented in cdk and not in ambit
    	return SubStructureSearchEngine( _smi, _sma, "cdk", _aromatic, _verbose );
    }

	/**
	 * Handle logical operators in a smarts class and uses atom by atom search (ABAS) to assign a compound.
	 * 
	 * @param _smiles
	 * @param smartsList
	 * @param module
	 * 
	 * @return
	 */
	public static boolean assignChemaxon( String _smiles, List<String> _smartsList, String _module, 
																boolean _aromatic, boolean _verbose ) {
		try {
			List<String> ORResponse = new ArrayList<String>();
			List<String> NOTResponse = new ArrayList<String>();
			
			boolean isORstructure = false;
			boolean isNOTstructure = false;
			String trueResponse = new String( "true" );
			int size = _smartsList.size();
		
			/*
			 * A class consists of smarts set.
			 * Smarts set is a collection of smarts query or smarts structures with logical operators
			 * the following additional logic has been added to for CDK smarts:
			 * 
			 * AND structures: query1XXXquery2, query1.query2. Both dots or XXX separated smarts are taken apart and we check each of them individually and only when all are ok then the whole expression is considered as AND connected structure.
			 * POLY structures: 3ZZZquery2 query 2 must occor 3 times
			 * NOT structures: !query, !query1XXXquery2, !query1.query2
			 * OR structures: queries without NOT symol, and AND structures after resolving XXX, WWW, dots.
			 */
			for ( int i=0; i < _smartsList.size(); i++ ) {
				
				String Qry = _smartsList.get(i);
				String Qry2 = Qry.replaceAll( "^!", "" );
				
	    		if ( Qry2.equals(Qry) ) {
	    			
	    			/*
	    			 * OR and AND structures. 
	    			 */
	    			String Qry1 = Qry.replaceAll( "XXX", "" );
	    			
	    			if( Qry1.equals( Qry ) ) {
	    				/*
	    				 * no patterns found, OR structure
	    				 */
		    			int rspd = SubStructureSearchEngine( _smiles, Qry, _module, _aromatic, _verbose );
		    			if ( rspd >0 ) ORResponse.add( "true" );
	    				else ORResponse.add( "false" );
	    				isORstructure = true;
		    			//System.out.println("sss: " + " OR " + rspd + " " + Qry);
	    		
	    			} else {
		    			/*
		    			 * XXX patterns found, AND structure
		    			 */
	    				List<String> newQryList = Segmenter( Qry, "XXX" );
	    				List<String> AndConnectedList= new ArrayList<String>();
	    				for ( int j = 0; j < newQryList.size(); j++ ) {
	    					String newQry = newQryList.get(j);
	    					int rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
	    					//System.out.println("sss: " + rspd + " " + newQry);
	    					if ( rspd >0 ) AndConnectedList.add("true");
	    					else AndConnectedList.add("false");
	    				}
	    				if ( IsAndConnected(AndConnectedList) ) { 
	    					String rspdS = AndConnectedList.get(0);
	    					ORResponse.add( rspdS );
	    					isORstructure = true;
	    					//System.out.println("AND with XXX " + rspd + " " + Qry);
	    				}
	    			}
	    		} else {
	    			/**
					 *  OR and AND structures preceeded by NOT
					 */
	    			String Qry1 = Qry2.replaceAll("XXX", "");
	    			
	    			if( Qry1.equals( Qry2 ) ) {
	    				/**
	    				 * no patterns found, OR structure preceeded by NOT
	    				 */
	    				int rspd1 = SubStructureSearchEngine( _smiles, Qry2, _module, _aromatic, _verbose );
	    				if ( rspd1 >0 ) NOTResponse.add("false");
    					else NOTResponse.add("true");
		    			isNOTstructure = true;
	    			} else {
		    			/**
		    			 * XXX patterns found, AND structure preceeded by NOT
		    			 */
	    				List<String> newQryList = Segmenter( Qry2, "XXX");
	    				List<String> AndConnectedList= new ArrayList<String>();
	    				
	    				for ( int j = 0; j<newQryList.size(); j++ ) {
	    					String newQry = newQryList.get(j);
	    					int rspd = SubStructureSearchEngine( _smiles, newQry, _module, _aromatic, _verbose );
	    					if ( rspd >0 ) AndConnectedList.add("true");
	    					else AndConnectedList.add("false");
	    				}
	    				if ( IsAndConnected(AndConnectedList ) ) { 
	    					String rspd1 = AndConnectedList.get(0);
	    					if ( rspd1.equals( "true" ) ) NOTResponse.add( "false" );
		             	    else NOTResponse.add( "true" );
	    					isNOTstructure = true;
	    				}
	    			}
	    		}
			}
			
			/*
			 *  a substructure match is considered, if atleast one OR structure is a match in the smarts set
			 * and all NOT structures in the smarts set did not match
			 */
			boolean asgn1 = IsOrConnected( ORResponse, "true" ) ;
			
			/*
			 * two scenarios in NOT structures: 1. only NOT structure 2. both NOT and OR structures
			 * scenario 1: only NOT
			 */
			if( isNOTstructure ) {
				boolean asgn2 = false;
				boolean check = IsAndConnected( NOTResponse );
				if ( check ) asgn2 = Boolean.valueOf ( NOTResponse.get(0) );
				
				if( !isORstructure ) {
					if( asgn2 ) return true;
					else return false;
				}
				
				/*
				 * two scenarios in NOT structures: 1. only NOT structure 2. both NOT and OR structures
				 * scenario 2: both NOT and OR
				 */
				if( isORstructure ) {
		    		if ( asgn2 ) {
		    			if ( asgn1 ) return true;
		    			else return false;
		    		}
		    	}
			}
			
			/*
			 * scenario 3: no NOT structures only OR structures
			 */
			else if ( isORstructure ) {
		    		//System.out.println("found OR structures");
		    		if ( asgn1 ) return true;
		    		else return false;
		    }
				
		} catch (Exception e) {
			System.out.println("assigner subroutine error: " + e);
		}
		return false;
	}
	
    /**
	 * the ABAS response for each smarts in the smarts set are AND connected or not
	 */
	public static boolean IsAndConnected( List<String> list ) {
		try {
		    for ( String s : list) {
		        if ( !s.equals( list.get(0) ) )
		            return false;
		    }
		    return true;
		} catch ( Exception e ) {
			System.out.println( "error IsAndConnected: " + e );
		}
		return false;
	}
	
    /**
	 * the ABAS response for each smarts in the smarts set are IS connected or not 
	 */
	public static boolean IsOrConnected(List<String> list, String toCheckValue) {
		try {
	        return list.contains( toCheckValue );
		} catch ( Exception e ) {
			System.out.println( "error IsOrConnected: " + e );
		}
		return false;
	}

    /**
	 * The AND connected smarts is represented as XXX or dots in the ontology. 
	 */
    public static List<String> Segmenter( String QueryWithPattern, String patternType ) {
		try {
			List<String> newQryList = Lists.newArrayList( Splitter.on( patternType ).trimResults().omitEmptyStrings().splitToList(QueryWithPattern));
			return newQryList;
		} catch( Exception e ) {
			System.out.println( "error smarts segmenter: " + e );
		}
		return null;
	}
    
	/**
	 * @param smi
	 * @param sma
	 * @param module
	 * @return true or false after processing POLY structures with keywords such as 1. EXACT, 2. MORE and 3. queries containing no keywords. 
	 * @throws EmptyMoleculeException
	 * @throws Exception
	 */
	public static int checkQueryMultiplicity( String _smi, String _sma, String _module, boolean _aromatic,  boolean _verbose ) throws Exception {
		try  {
			// multiplicity keywords are only used for Cdk or Ambit handling.
			String exact = _sma.replaceAll( "EXACT", "" );
			String more = _sma.replaceAll( "MORE", "" );
			
			if ( exact.equals( _sma ) ) {
				//does not contain EXACT keyword
				if ( more.equals( _sma ) ) {
					//does not contain MORE keyword
	 				return StructureSearchEngine.searchBySubstructure( _smi, _sma, _module, _aromatic, _verbose );
				} else {
					//contains keyword MORE
					List<String> listSmarts = Segmenter( _sma, "MORE" );
					String query = listSmarts.get(1);
					String threshold = listSmarts.get(0);
					int cnt = StructureSearchEngine.searchBySubstructureAmbitAllInstances( _smi, query, _verbose );
					if (cnt >= Integer.parseInt( threshold ) ) return 1;
					else  return 0;
				}
	 		} else {
				//contains keyword EXACT
				List<String> listSmarts = Segmenter( _sma, "EXACT" );
				String query = listSmarts.get(1);
				String threshold = listSmarts.get(0);
				int cnt = StructureSearchEngine.searchBySubstructureAmbitAllInstances( _smi, query, _verbose );
				if ( cnt == Integer.parseInt(threshold) ) return 1;
				else return 0;
			}
		} catch ( Exception e ) {
			LOG.info( "ERROR: error in processing smarts multiplicity: " + e ) ;
			return -1;
		}
	}
	
	/**
	 * perform atom-by-atom-search (ABAS) given a query and a target 
	 */
	private static int SubStructureSearchEngine( String _target, String _query, String _module, 
																		boolean _aromatic, boolean _verbose ) {
		try {
			return StructureSearchEngine.searchBySubstructure( _target, _query, _module, _aromatic, _verbose ) ;
		} catch (Exception e)  {
			LOG.info( "ERROR: SubStructureSearchEngine Error "+e);
		}
		return -1;
	}
    
}
