package COBRApyBasedFBA::COBRApyBasedFBAClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

COBRApyBasedFBA::COBRApyBasedFBAClient

=head1 DESCRIPTION


A KBase module: COBRApyBasedFBA


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => COBRApyBasedFBA::COBRApyBasedFBAClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_fba_pipeline

  $results = $obj->run_fba_pipeline($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a COBRApyBasedFBA.RunFBAPipelineParams
$results is a COBRApyBasedFBA.RunFBAPipelineResults
RunFBAPipelineParams is a reference to a hash where the following keys are defined:
	fbamodel_id has a value which is a COBRApyBasedFBA.fbamodel_id
	fbamodel_workspace has a value which is a COBRApyBasedFBA.workspace_name
	media_id has a value which is a COBRApyBasedFBA.media_id
	media_workspace has a value which is a COBRApyBasedFBA.workspace_name
	target_reaction has a value which is a COBRApyBasedFBA.reaction_id
	fba_output_id has a value which is a COBRApyBasedFBA.fba_id
	workspace has a value which is a COBRApyBasedFBA.workspace_name
	fva has a value which is a COBRApyBasedFBA.bool
	minimize_flux has a value which is a COBRApyBasedFBA.bool
	simulate_ko has a value which is a COBRApyBasedFBA.bool
	all_reversible has a value which is a COBRApyBasedFBA.bool
	feature_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.feature_id
	reaction_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.reaction_id
	media_supplement_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.compound_id
	objective_fraction has a value which is a float
	max_c_uptake has a value which is a float
	max_n_uptake has a value which is a float
	max_p_uptake has a value which is a float
	max_s_uptake has a value which is a float
	max_o_uptake has a value which is a float
	default_max_uptake has a value which is a float
fbamodel_id is a string
workspace_name is a string
media_id is a string
reaction_id is a string
fba_id is a string
bool is an int
feature_id is a string
compound_id is a string
RunFBAPipelineResults is a reference to a hash where the following keys are defined:
	new_fba_ref has a value which is a COBRApyBasedFBA.ws_fba_id
	objective has a value which is an int
	report_name has a value which is a string
	report_ref has a value which is a COBRApyBasedFBA.ws_report_id
ws_fba_id is a string
ws_report_id is a string

</pre>

=end html

=begin text

$params is a COBRApyBasedFBA.RunFBAPipelineParams
$results is a COBRApyBasedFBA.RunFBAPipelineResults
RunFBAPipelineParams is a reference to a hash where the following keys are defined:
	fbamodel_id has a value which is a COBRApyBasedFBA.fbamodel_id
	fbamodel_workspace has a value which is a COBRApyBasedFBA.workspace_name
	media_id has a value which is a COBRApyBasedFBA.media_id
	media_workspace has a value which is a COBRApyBasedFBA.workspace_name
	target_reaction has a value which is a COBRApyBasedFBA.reaction_id
	fba_output_id has a value which is a COBRApyBasedFBA.fba_id
	workspace has a value which is a COBRApyBasedFBA.workspace_name
	fva has a value which is a COBRApyBasedFBA.bool
	minimize_flux has a value which is a COBRApyBasedFBA.bool
	simulate_ko has a value which is a COBRApyBasedFBA.bool
	all_reversible has a value which is a COBRApyBasedFBA.bool
	feature_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.feature_id
	reaction_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.reaction_id
	media_supplement_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.compound_id
	objective_fraction has a value which is a float
	max_c_uptake has a value which is a float
	max_n_uptake has a value which is a float
	max_p_uptake has a value which is a float
	max_s_uptake has a value which is a float
	max_o_uptake has a value which is a float
	default_max_uptake has a value which is a float
fbamodel_id is a string
workspace_name is a string
media_id is a string
reaction_id is a string
fba_id is a string
bool is an int
feature_id is a string
compound_id is a string
RunFBAPipelineResults is a reference to a hash where the following keys are defined:
	new_fba_ref has a value which is a COBRApyBasedFBA.ws_fba_id
	objective has a value which is an int
	report_name has a value which is a string
	report_ref has a value which is a COBRApyBasedFBA.ws_report_id
ws_fba_id is a string
ws_report_id is a string


=end text

=item Description

Run flux balance analysis and return ID of FBA object with results

=back

=cut

 sub run_fba_pipeline
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_fba_pipeline (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_fba_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_fba_pipeline');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "COBRApyBasedFBA.run_fba_pipeline",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_fba_pipeline',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_fba_pipeline",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_fba_pipeline',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "COBRApyBasedFBA.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "COBRApyBasedFBA.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_fba_pipeline',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_fba_pipeline",
            status_line => $self->{client}->status_line,
            method_name => 'run_fba_pipeline',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for COBRApyBasedFBA::COBRApyBasedFBAClient\n";
    }
    if ($sMajor == 0) {
        warn "COBRApyBasedFBA::COBRApyBasedFBAClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 bool

=over 4



=item Description

A binary boolean


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 media_id

=over 4



=item Description

A string representing a Media id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 fbamodel_id

=over 4



=item Description

A string representing a FBAModel id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 fba_id

=over 4



=item Description

A string representing a FBA id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 reaction_id

=over 4



=item Description

A string representing a reaction id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_id

=over 4



=item Description

A string representing a feature id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 compound_id

=over 4



=item Description

A string representing a compound id.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 workspace_name

=over 4



=item Description

A string representing a workspace name.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_fbamodel_id

=over 4



=item Description

The workspace ID for a FBAModel data object.
@id ws KBaseFBA.FBAModel


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_fba_id

=over 4



=item Description

The workspace ID for a FBA data object.
@id ws KBaseFBA.FBA


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 ws_report_id

=over 4



=item Description

The workspace ID for a Report object
@id ws KBaseReport.Report


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 RunFBAPipelineParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
fbamodel_id has a value which is a COBRApyBasedFBA.fbamodel_id
fbamodel_workspace has a value which is a COBRApyBasedFBA.workspace_name
media_id has a value which is a COBRApyBasedFBA.media_id
media_workspace has a value which is a COBRApyBasedFBA.workspace_name
target_reaction has a value which is a COBRApyBasedFBA.reaction_id
fba_output_id has a value which is a COBRApyBasedFBA.fba_id
workspace has a value which is a COBRApyBasedFBA.workspace_name
fva has a value which is a COBRApyBasedFBA.bool
minimize_flux has a value which is a COBRApyBasedFBA.bool
simulate_ko has a value which is a COBRApyBasedFBA.bool
all_reversible has a value which is a COBRApyBasedFBA.bool
feature_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.feature_id
reaction_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.reaction_id
media_supplement_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.compound_id
objective_fraction has a value which is a float
max_c_uptake has a value which is a float
max_n_uptake has a value which is a float
max_p_uptake has a value which is a float
max_s_uptake has a value which is a float
max_o_uptake has a value which is a float
default_max_uptake has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
fbamodel_id has a value which is a COBRApyBasedFBA.fbamodel_id
fbamodel_workspace has a value which is a COBRApyBasedFBA.workspace_name
media_id has a value which is a COBRApyBasedFBA.media_id
media_workspace has a value which is a COBRApyBasedFBA.workspace_name
target_reaction has a value which is a COBRApyBasedFBA.reaction_id
fba_output_id has a value which is a COBRApyBasedFBA.fba_id
workspace has a value which is a COBRApyBasedFBA.workspace_name
fva has a value which is a COBRApyBasedFBA.bool
minimize_flux has a value which is a COBRApyBasedFBA.bool
simulate_ko has a value which is a COBRApyBasedFBA.bool
all_reversible has a value which is a COBRApyBasedFBA.bool
feature_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.feature_id
reaction_ko_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.reaction_id
media_supplement_list has a value which is a reference to a list where each element is a COBRApyBasedFBA.compound_id
objective_fraction has a value which is a float
max_c_uptake has a value which is a float
max_n_uptake has a value which is a float
max_p_uptake has a value which is a float
max_s_uptake has a value which is a float
max_o_uptake has a value which is a float
default_max_uptake has a value which is a float


=end text

=back



=head2 RunFBAPipelineResults

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
new_fba_ref has a value which is a COBRApyBasedFBA.ws_fba_id
objective has a value which is an int
report_name has a value which is a string
report_ref has a value which is a COBRApyBasedFBA.ws_report_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
new_fba_ref has a value which is a COBRApyBasedFBA.ws_fba_id
objective has a value which is an int
report_name has a value which is a string
report_ref has a value which is a COBRApyBasedFBA.ws_report_id


=end text

=back



=cut

package COBRApyBasedFBA::COBRApyBasedFBAClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
