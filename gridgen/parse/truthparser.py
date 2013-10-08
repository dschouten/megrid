from ROOT import TLorentzVector as tlv

# ------------------------------------------------------------------------------

PDG = { 'W' : 24, 'W-' : 24, 'W+' : -24,
        'e' : 11, 'e-' : 11, 'e+' : -11,
        'mu' : 13, 'mu-' : 13, 'mu+' : -13,
        'tau' : 15, 'tau-' : 15, 'tau+' : -15,
        've' : 12, '~ve' : -12,
        'vmu' : 14, '~vmu' : -14,
        'vtau' : 16, '~vtau' : -16,
        'd' : 1,  'u': 2,  's': 3,  'c': 4,  'b': 5,  't': 6,
        '~d':-1, '~u':-2, '~s':-3, '~c':-4, '~b':-5, '~t':-6 }

tlv_lp   = tlv()
tlv_vl   = tlv()
tlv_lm   = tlv()
tlv_vr   = tlv()
tlv_wp   = tlv()
tlv_wm   = tlv()
tlv_b    = tlv()
tlv_bbar = tlv()
tlv_t    = tlv()
tlv_tbar = tlv()

# ------------------------------------------------------------------------------

def parse_ww_pythia( itree, lepdecays=( PDG['e'],PDG['mu'],PDG['tau'] ), version='PYTHIA8', indices=[], status=20 ):
    global tlv_wm, tlv_lm, tlv_vr, tlv_wp, tlv_lp, tlv_vl
    ilp = -1
    ilm = -1
    ivl = -1
    ivr = -1
    iwp = -1
    iwm = -1
    if len( indices ) == 0:
        indices = xrange( itree.mc_n )
    iversion=int(version[-1])
    for imc in indices:        
        if ilp<0 or ilm<0 or ivr<0 or ivl<0 or iwp<0 or iwm<0:
            if abs(itree.mc_pdgId[imc]) == PDG['W'] and ( (iversion==8 and itree.mc_status[imc] > status) or
                                                          (iversion==6 and itree.mc_barcode[imc] < 20) ):
                haslepdecay = 0
                for ichild in itree.mc_child_index[imc]:
                    if ilp<0 or ilm<0:
                        if abs(itree.mc_pdgId[ichild]) in lepdecays:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ilm = ichild
                                else:                       ilp = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in lepdecays:
                                        if itree.mc_pdgId[jmc] > 0: ilm = jmc
                                        else:                       ilp = jmc
                                        jmc = kmc
                    if ivr<0 or ivl<0:
                        if abs(itree.mc_pdgId[ichild]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ivl = ichild
                                else:                       ivr = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                                        if itree.mc_pdgId[jmc] > 0: ivl = jmc
                                        else:                       ivr = jmc
                                        jmc = kmc
                if haslepdecay > 1:
                    if itree.mc_pdgId[imc] > 0:
                        iwm = imc
                    else:
                        iwp = imc
                            
        if ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0:
            break

    if not( ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0 ):
        return None, None, None, None, None, None, -1, -1

    tlv_lp.SetPtEtaPhiE( itree.mc_pt[ilp], itree.mc_eta[ilp], itree.mc_phi[ilp], itree.mc_E[ilp] )
    tlv_vl.SetPtEtaPhiE( itree.mc_pt[ivl], itree.mc_eta[ivl], itree.mc_phi[ivl], itree.mc_E[ivl] )
    tlv_lm.SetPtEtaPhiE( itree.mc_pt[ilm], itree.mc_eta[ilm], itree.mc_phi[ilm], itree.mc_E[ilm] )
    tlv_vr.SetPtEtaPhiE( itree.mc_pt[ivr], itree.mc_eta[ivr], itree.mc_phi[ivr], itree.mc_E[ivr] )

    tlv_wp.SetPtEtaPhiE( itree.mc_pt[iwp], itree.mc_eta[iwp], itree.mc_phi[iwp], itree.mc_E[iwp] )
    tlv_wm.SetPtEtaPhiE( itree.mc_pt[iwm], itree.mc_eta[iwm], itree.mc_phi[iwm], itree.mc_E[iwm] )

    return (tlv_wm, tlv_lm, tlv_vr, tlv_wp, tlv_lp, tlv_vl, itree.mc_pdgId[ilm], itree.mc_pdgId[ilp] )

# ------------------------------------------------------------------------------

parse_ww_herwig   = None
parse_ww_herwigpp = None

# ------------------------------------------------------------------------------

def parse_ttbar( itree, lepdecays=( PDG['e'],PDG['mu'],PDG['tau'] ), version=None, indices=[], status=20 ):
    global tlv_b, tlv_bbar, tlv_wm, tlv_lm, tlv_vr, tlv_wp, tlv_lp, tlv_vl
    ilp = -1
    ilm = -1
    ivl = -1
    ivr = -1
    iwp = -1
    iwm = -1

    it    = -1
    itbar = -1
    ib    = -1
    ibbar = -1

    if len(indices)==0:
        indices = xrange( itree.mc_n )
    
    # if len( indices ) != 0: t_indices = indices[:]
    # if version != None:
    #     iversion = int(version[-1])
    
    for imc in indices:        
        if it < 0 or itbar < 0 or ib < 0 or ibbar < 0:
            #
            # find W->t,b vertex, define W vertices
            #
            if abs(itree.mc_pdgId[imc]) == PDG['t'] and itree.mc_child_index[imc].size() != 0:
                hasbdecay = 0
                haswdecay = 0
                for jmc in itree.mc_child_index[imc]:
                    if abs(itree.mc_pdgId[jmc]) == PDG['b']:
                        hasbdecay += 1
                    if abs(itree.mc_pdgId[jmc]) == PDG['W']:
                        haswdecay += 1
                if hasbdecay==1 and haswdecay==1:
                    for jmc in itree.mc_child_index[imc]:
                        if itree.mc_pdgId[jmc] == PDG['b']:
                            tlv_b.SetPtEtaPhiM( itree.mc_pt[jmc], itree.mc_eta[jmc], itree.mc_phi[jmc], itree.mc_m[jmc] )
                        if itree.mc_pdgId[jmc] == PDG['~b']:
                            tlv_bbar.SetPtEtaPhiM( itree.mc_pt[jmc], itree.mc_eta[jmc], itree.mc_phi[jmc], itree.mc_m[jmc] )
                        if itree.mc_pdgId[jmc] == PDG['W-']:
                            ## tlv_wm.SetPtEtaPhiM( itree.mc_pt[jmc], itree.mc_eta[jmc], itree.mc_phi[jmc], itree.mc_m[jmc] )
                            ## w_indices += [jmc]
                            pass
                        if itree.mc_pdgId[jmc] == PDG['W+']:
                            ## tlv_wp.SetPtEtaPhiM( itree.mc_pt[jmc], itree.mc_eta[jmc], itree.mc_phi[jmc], itree.mc_m[jmc] )
                            ## w_indices += [jmc]
                            pass
    
    for imc in indices:
        if ilp<0 or ilm<0 or ivr<0 or ivl<0 or iwp<0 or iwm<0:
            if abs(itree.mc_pdgId[imc]) == PDG['W']:
                #
                # loop over W decays ...
                #
                haslepdecay = 0
                for ichild in itree.mc_child_index[imc]:
                    if ilp<0 or ilm<0:
                        if abs(itree.mc_pdgId[ichild]) in lepdecays:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ilm = ichild
                                else:                       ilp = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in lepdecays:
                                        if itree.mc_pdgId[jmc] > 0: ilm = jmc
                                        else:                       ilp = jmc
                                        jmc = kmc
                    if ivr<0 or ivl<0:
                        if abs(itree.mc_pdgId[ichild]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ivl = ichild
                                else:                       ivr = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                                        if itree.mc_pdgId[jmc] > 0: ivl = jmc
                                        else:                       ivr = jmc
                                        jmc = kmc
                if haslepdecay > 1:
                    if itree.mc_pdgId[imc] > 0:
                        iwm = imc
                    else:
                        iwp = imc
                            
        if ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0:
            break

    if not( ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0 ):
        return None, None, None, None, None, None, -1, -1

    tlv_lp.SetPtEtaPhiE( itree.mc_pt[ilp], itree.mc_eta[ilp], itree.mc_phi[ilp], itree.mc_E[ilp] )
    tlv_vl.SetPtEtaPhiE( itree.mc_pt[ivl], itree.mc_eta[ivl], itree.mc_phi[ivl], itree.mc_E[ivl] )
    tlv_lm.SetPtEtaPhiE( itree.mc_pt[ilm], itree.mc_eta[ilm], itree.mc_phi[ilm], itree.mc_E[ilm] )
    tlv_vr.SetPtEtaPhiE( itree.mc_pt[ivr], itree.mc_eta[ivr], itree.mc_phi[ivr], itree.mc_E[ivr] )

    return (tlv_b, tlv_lm, tlv_vr, tlv_bbar, tlv_lp, tlv_vl, itree.mc_pdgId[ilm], itree.mc_pdgId[ilp])

# ------------------------------------------------------------------------------

parse_ttbar_herwig = parse_ttbar
parse_ttbar_pythia = parse_ttbar

# ------------------------------------------------------------------------------

def parse_wt( itree, lepdecays=( PDG['e'],PDG['mu'],PDG['tau'] ), version=None, indices=[], status=20 ):
    global tlv_b, tlv_t, tlv_wm, tlv_lm, tlv_vr, tlv_wp, tlv_lp, tlv_vl
    ilp = -1
    ilm = -1
    ivl = -1
    ivr = -1
    iwp = -1
    iwm = -1
    it  = -1
    ib  = -1

    if len(indices)==0:
        indices = xrange( itree.mc_n )
    
    # if len( indices ) != 0: t_indices = indices[:]
    # if version != None:
    #     iversion = int(version[-1])
    
    for imc in indices:        
        if it < 0 or ib < 0:
            #
            # find W->t,b vertex
            #
            if abs(itree.mc_pdgId[imc]) == PDG['t'] and itree.mc_child_index[imc].size() != 0:
                hasbdecay = 0
                haswdecay = 0
                for jmc in itree.mc_child_index[imc]:
                    if abs(itree.mc_pdgId[jmc]) == PDG['b']:
                        hasbdecay += 1
                    if abs(itree.mc_pdgId[jmc]) == PDG['W']:
                        haswdecay += 1
                if hasbdecay==1 and haswdecay==1:
                    it = imc
                    tlv_t.SetPtEtaPhiM( itree.mc_pt[imc], itree.mc_eta[imc], itree.mc_phi[imc], itree.mc_m[imc] )
                    for jmc in itree.mc_child_index[imc]:
                        if abs(itree.mc_pdgId[jmc]) == PDG['b']:
                            ib = jmc
                            tlv_b.SetPtEtaPhiM( itree.mc_pt[jmc], itree.mc_eta[jmc], itree.mc_phi[jmc], itree.mc_m[jmc] )
    
    for imc in indices:
        if ilp<0 or ilm<0 or ivr<0 or ivl<0 or iwp<0 or iwm<0:
            if abs(itree.mc_pdgId[imc]) == PDG['W']:
                #
                # loop over W decays ...
                #
                haslepdecay = 0
                for ichild in itree.mc_child_index[imc]:
                    if ilp<0 or ilm<0:
                        if abs(itree.mc_pdgId[ichild]) in lepdecays:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ilm = ichild
                                else:                       ilp = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in lepdecays:
                                        if itree.mc_pdgId[jmc] > 0: ilm = jmc
                                        else:                       ilp = jmc
                                        jmc = kmc
                    if ivr<0 or ivl<0:
                        if abs(itree.mc_pdgId[ichild]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                            haslepdecay += 1
                            decayed = False
                            jmc = ichild
                            while not decayed:
                                if itree.mc_pdgId[jmc] > 0: ivl = ichild
                                else:                       ivr = ichild
                                decayed = True
                                for kmc in itree.mc_child_index[jmc]:
                                    decayed = False
                                    if abs(itree.mc_pdgId[jmc]) in [PDG['ve'],PDG['vmu'],PDG['vtau']]:
                                        if itree.mc_pdgId[jmc] > 0: ivl = jmc
                                        else:                       ivr = jmc
                                        jmc = kmc
                if haslepdecay > 1:
                    if itree.mc_pdgId[imc] > 0:
                        iwm = imc
                    else:
                        iwp = imc
                            
        if ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0:
            break

    if not( ilp >= 0 and ivr >= 0 and ilm >= 0 and ivl >= 0 and iwm >= 0 and iwp >= 0 ):
        return None, None, None, None, None, None, -1, -1

    tlv_lp.SetPtEtaPhiE( itree.mc_pt[ilp], itree.mc_eta[ilp], itree.mc_phi[ilp], itree.mc_E[ilp] )
    tlv_vl.SetPtEtaPhiE( itree.mc_pt[ivl], itree.mc_eta[ivl], itree.mc_phi[ivl], itree.mc_E[ivl] )
    tlv_lm.SetPtEtaPhiE( itree.mc_pt[ilm], itree.mc_eta[ilm], itree.mc_phi[ilm], itree.mc_E[ilm] )
    tlv_vr.SetPtEtaPhiE( itree.mc_pt[ivr], itree.mc_eta[ivr], itree.mc_phi[ivr], itree.mc_E[ivr] )

    return (tlv_b, tlv_lm, tlv_vr, tlv_lp, tlv_vl, itree.mc_pdgId[ilm], itree.mc_pdgId[ilp], itree.mc_pdgId[ib])

# ------------------------------------------------------------------------------

parse_wt_herwig = parse_wt
parse_wt_pythia = parse_wt
