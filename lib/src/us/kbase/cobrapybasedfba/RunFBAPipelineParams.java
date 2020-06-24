
package us.kbase.cobrapybasedfba;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: RunFBAPipelineParams</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "fbamodel_id",
    "fbamodel_workspace",
    "media_id",
    "media_workspace",
    "target_reaction",
    "fba_output_id",
    "workspace",
    "fva",
    "minimize_flux",
    "simulate_ko",
    "all_reversible",
    "feature_ko_list",
    "reaction_ko_list",
    "media_supplement_list",
    "objective_fraction",
    "max_c_uptake",
    "max_n_uptake",
    "max_p_uptake",
    "max_s_uptake",
    "max_o_uptake",
    "default_max_uptake"
})
public class RunFBAPipelineParams {

    @JsonProperty("fbamodel_id")
    private java.lang.String fbamodelId;
    @JsonProperty("fbamodel_workspace")
    private java.lang.String fbamodelWorkspace;
    @JsonProperty("media_id")
    private java.lang.String mediaId;
    @JsonProperty("media_workspace")
    private java.lang.String mediaWorkspace;
    @JsonProperty("target_reaction")
    private java.lang.String targetReaction;
    @JsonProperty("fba_output_id")
    private java.lang.String fbaOutputId;
    @JsonProperty("workspace")
    private java.lang.String workspace;
    @JsonProperty("fva")
    private Long fva;
    @JsonProperty("minimize_flux")
    private Long minimizeFlux;
    @JsonProperty("simulate_ko")
    private Long simulateKo;
    @JsonProperty("all_reversible")
    private Long allReversible;
    @JsonProperty("feature_ko_list")
    private List<String> featureKoList;
    @JsonProperty("reaction_ko_list")
    private List<String> reactionKoList;
    @JsonProperty("media_supplement_list")
    private List<String> mediaSupplementList;
    @JsonProperty("objective_fraction")
    private Double objectiveFraction;
    @JsonProperty("max_c_uptake")
    private Double maxCUptake;
    @JsonProperty("max_n_uptake")
    private Double maxNUptake;
    @JsonProperty("max_p_uptake")
    private Double maxPUptake;
    @JsonProperty("max_s_uptake")
    private Double maxSUptake;
    @JsonProperty("max_o_uptake")
    private Double maxOUptake;
    @JsonProperty("default_max_uptake")
    private Double defaultMaxUptake;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("fbamodel_id")
    public java.lang.String getFbamodelId() {
        return fbamodelId;
    }

    @JsonProperty("fbamodel_id")
    public void setFbamodelId(java.lang.String fbamodelId) {
        this.fbamodelId = fbamodelId;
    }

    public RunFBAPipelineParams withFbamodelId(java.lang.String fbamodelId) {
        this.fbamodelId = fbamodelId;
        return this;
    }

    @JsonProperty("fbamodel_workspace")
    public java.lang.String getFbamodelWorkspace() {
        return fbamodelWorkspace;
    }

    @JsonProperty("fbamodel_workspace")
    public void setFbamodelWorkspace(java.lang.String fbamodelWorkspace) {
        this.fbamodelWorkspace = fbamodelWorkspace;
    }

    public RunFBAPipelineParams withFbamodelWorkspace(java.lang.String fbamodelWorkspace) {
        this.fbamodelWorkspace = fbamodelWorkspace;
        return this;
    }

    @JsonProperty("media_id")
    public java.lang.String getMediaId() {
        return mediaId;
    }

    @JsonProperty("media_id")
    public void setMediaId(java.lang.String mediaId) {
        this.mediaId = mediaId;
    }

    public RunFBAPipelineParams withMediaId(java.lang.String mediaId) {
        this.mediaId = mediaId;
        return this;
    }

    @JsonProperty("media_workspace")
    public java.lang.String getMediaWorkspace() {
        return mediaWorkspace;
    }

    @JsonProperty("media_workspace")
    public void setMediaWorkspace(java.lang.String mediaWorkspace) {
        this.mediaWorkspace = mediaWorkspace;
    }

    public RunFBAPipelineParams withMediaWorkspace(java.lang.String mediaWorkspace) {
        this.mediaWorkspace = mediaWorkspace;
        return this;
    }

    @JsonProperty("target_reaction")
    public java.lang.String getTargetReaction() {
        return targetReaction;
    }

    @JsonProperty("target_reaction")
    public void setTargetReaction(java.lang.String targetReaction) {
        this.targetReaction = targetReaction;
    }

    public RunFBAPipelineParams withTargetReaction(java.lang.String targetReaction) {
        this.targetReaction = targetReaction;
        return this;
    }

    @JsonProperty("fba_output_id")
    public java.lang.String getFbaOutputId() {
        return fbaOutputId;
    }

    @JsonProperty("fba_output_id")
    public void setFbaOutputId(java.lang.String fbaOutputId) {
        this.fbaOutputId = fbaOutputId;
    }

    public RunFBAPipelineParams withFbaOutputId(java.lang.String fbaOutputId) {
        this.fbaOutputId = fbaOutputId;
        return this;
    }

    @JsonProperty("workspace")
    public java.lang.String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(java.lang.String workspace) {
        this.workspace = workspace;
    }

    public RunFBAPipelineParams withWorkspace(java.lang.String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonProperty("fva")
    public Long getFva() {
        return fva;
    }

    @JsonProperty("fva")
    public void setFva(Long fva) {
        this.fva = fva;
    }

    public RunFBAPipelineParams withFva(Long fva) {
        this.fva = fva;
        return this;
    }

    @JsonProperty("minimize_flux")
    public Long getMinimizeFlux() {
        return minimizeFlux;
    }

    @JsonProperty("minimize_flux")
    public void setMinimizeFlux(Long minimizeFlux) {
        this.minimizeFlux = minimizeFlux;
    }

    public RunFBAPipelineParams withMinimizeFlux(Long minimizeFlux) {
        this.minimizeFlux = minimizeFlux;
        return this;
    }

    @JsonProperty("simulate_ko")
    public Long getSimulateKo() {
        return simulateKo;
    }

    @JsonProperty("simulate_ko")
    public void setSimulateKo(Long simulateKo) {
        this.simulateKo = simulateKo;
    }

    public RunFBAPipelineParams withSimulateKo(Long simulateKo) {
        this.simulateKo = simulateKo;
        return this;
    }

    @JsonProperty("all_reversible")
    public Long getAllReversible() {
        return allReversible;
    }

    @JsonProperty("all_reversible")
    public void setAllReversible(Long allReversible) {
        this.allReversible = allReversible;
    }

    public RunFBAPipelineParams withAllReversible(Long allReversible) {
        this.allReversible = allReversible;
        return this;
    }

    @JsonProperty("feature_ko_list")
    public List<String> getFeatureKoList() {
        return featureKoList;
    }

    @JsonProperty("feature_ko_list")
    public void setFeatureKoList(List<String> featureKoList) {
        this.featureKoList = featureKoList;
    }

    public RunFBAPipelineParams withFeatureKoList(List<String> featureKoList) {
        this.featureKoList = featureKoList;
        return this;
    }

    @JsonProperty("reaction_ko_list")
    public List<String> getReactionKoList() {
        return reactionKoList;
    }

    @JsonProperty("reaction_ko_list")
    public void setReactionKoList(List<String> reactionKoList) {
        this.reactionKoList = reactionKoList;
    }

    public RunFBAPipelineParams withReactionKoList(List<String> reactionKoList) {
        this.reactionKoList = reactionKoList;
        return this;
    }

    @JsonProperty("media_supplement_list")
    public List<String> getMediaSupplementList() {
        return mediaSupplementList;
    }

    @JsonProperty("media_supplement_list")
    public void setMediaSupplementList(List<String> mediaSupplementList) {
        this.mediaSupplementList = mediaSupplementList;
    }

    public RunFBAPipelineParams withMediaSupplementList(List<String> mediaSupplementList) {
        this.mediaSupplementList = mediaSupplementList;
        return this;
    }

    @JsonProperty("objective_fraction")
    public Double getObjectiveFraction() {
        return objectiveFraction;
    }

    @JsonProperty("objective_fraction")
    public void setObjectiveFraction(Double objectiveFraction) {
        this.objectiveFraction = objectiveFraction;
    }

    public RunFBAPipelineParams withObjectiveFraction(Double objectiveFraction) {
        this.objectiveFraction = objectiveFraction;
        return this;
    }

    @JsonProperty("max_c_uptake")
    public Double getMaxCUptake() {
        return maxCUptake;
    }

    @JsonProperty("max_c_uptake")
    public void setMaxCUptake(Double maxCUptake) {
        this.maxCUptake = maxCUptake;
    }

    public RunFBAPipelineParams withMaxCUptake(Double maxCUptake) {
        this.maxCUptake = maxCUptake;
        return this;
    }

    @JsonProperty("max_n_uptake")
    public Double getMaxNUptake() {
        return maxNUptake;
    }

    @JsonProperty("max_n_uptake")
    public void setMaxNUptake(Double maxNUptake) {
        this.maxNUptake = maxNUptake;
    }

    public RunFBAPipelineParams withMaxNUptake(Double maxNUptake) {
        this.maxNUptake = maxNUptake;
        return this;
    }

    @JsonProperty("max_p_uptake")
    public Double getMaxPUptake() {
        return maxPUptake;
    }

    @JsonProperty("max_p_uptake")
    public void setMaxPUptake(Double maxPUptake) {
        this.maxPUptake = maxPUptake;
    }

    public RunFBAPipelineParams withMaxPUptake(Double maxPUptake) {
        this.maxPUptake = maxPUptake;
        return this;
    }

    @JsonProperty("max_s_uptake")
    public Double getMaxSUptake() {
        return maxSUptake;
    }

    @JsonProperty("max_s_uptake")
    public void setMaxSUptake(Double maxSUptake) {
        this.maxSUptake = maxSUptake;
    }

    public RunFBAPipelineParams withMaxSUptake(Double maxSUptake) {
        this.maxSUptake = maxSUptake;
        return this;
    }

    @JsonProperty("max_o_uptake")
    public Double getMaxOUptake() {
        return maxOUptake;
    }

    @JsonProperty("max_o_uptake")
    public void setMaxOUptake(Double maxOUptake) {
        this.maxOUptake = maxOUptake;
    }

    public RunFBAPipelineParams withMaxOUptake(Double maxOUptake) {
        this.maxOUptake = maxOUptake;
        return this;
    }

    @JsonProperty("default_max_uptake")
    public Double getDefaultMaxUptake() {
        return defaultMaxUptake;
    }

    @JsonProperty("default_max_uptake")
    public void setDefaultMaxUptake(Double defaultMaxUptake) {
        this.defaultMaxUptake = defaultMaxUptake;
    }

    public RunFBAPipelineParams withDefaultMaxUptake(Double defaultMaxUptake) {
        this.defaultMaxUptake = defaultMaxUptake;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((((((((((((((((((((((((((((((((((((((((((("RunFBAPipelineParams"+" [fbamodelId=")+ fbamodelId)+", fbamodelWorkspace=")+ fbamodelWorkspace)+", mediaId=")+ mediaId)+", mediaWorkspace=")+ mediaWorkspace)+", targetReaction=")+ targetReaction)+", fbaOutputId=")+ fbaOutputId)+", workspace=")+ workspace)+", fva=")+ fva)+", minimizeFlux=")+ minimizeFlux)+", simulateKo=")+ simulateKo)+", allReversible=")+ allReversible)+", featureKoList=")+ featureKoList)+", reactionKoList=")+ reactionKoList)+", mediaSupplementList=")+ mediaSupplementList)+", objectiveFraction=")+ objectiveFraction)+", maxCUptake=")+ maxCUptake)+", maxNUptake=")+ maxNUptake)+", maxPUptake=")+ maxPUptake)+", maxSUptake=")+ maxSUptake)+", maxOUptake=")+ maxOUptake)+", defaultMaxUptake=")+ defaultMaxUptake)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
