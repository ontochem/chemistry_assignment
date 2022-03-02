/**
 * Copyright OntoChem GmbH
 */
package com.ontochem.assignment;

import java.util.Locale;

/**
 * Constants for chemical libraries.
 *
 *
 * <h3>Changelog</h3>
 * <ul>
 *   <li>2022-02-25
 *     <ul>
 *       <li>initial version</li>
 *     </ul>
 *   </li>
 * </ul>
 */
public class ChemLib {

  public final static String CHEMLIB_AMBIT = "ambit";
  public final static String CHEMLIB_CDK   = "cdk";
  public final static String CHEMLIB_CA    = "chemaxon";

  private final static String[] KNOWN_CHEMLIBS = new String[] {
      CHEMLIB_AMBIT,
      CHEMLIB_CA,
      CHEMLIB_CDK
  };
  
  /**
   * Returns normalized chemical library name for provided name.
   * 
   * @param _name
   * 
   * @return  known normalized library name or <code>null</code>
   */
  public final static String resolveChemLib( String _name ) {
    if ( _name != null ) {
      final String nameLc = _name.toLowerCase( Locale.ENGLISH );
      for( String chemLib : KNOWN_CHEMLIBS ) {
        if ( chemLib.equals( nameLc ) ) {
          return chemLib;
        }
      }
    }
    return null;
  }
  
  /**
   * Returns <code>true</code> if provided name is a known chemical
   * library ({@link #resolveChemLib(String)} returns a non-null value).
   * 
   * @param _name
   * 
   * @return
   */
  public final static boolean isKnownChemLib( String _name ) {
    return resolveChemLib( _name ) != null;
  }
}
