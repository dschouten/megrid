void nevents( const char* F, const char* T )
{
  Int_t nEventsTotal;
  TFile* fp = TFile::Open( F, "read" );
  if( fp != NULL && !(fp->IsZombie()) )
  {
    TTree* tp = (TTree*)fp->Get( T );
    if( tp != NULL )
    {
      nEventsTotal = tp->GetEntries();
    }
  }
  std::cout << nEventsTotal << std::endl;
}
